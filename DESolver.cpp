//////////////////////////////////////////////////////////////////////
//
// Differential Evolution Test Program (MINIMIZATION)
// Based on algorithms developed by Dr. Rainer Storn & Kenneth Price
//
// Parallel version By:
//    Joao C. C. Henriques
//    joaochenriques@ist.utl.pt, joaochenriques@gmail.com
//    Institute of Mechanical Engineering - IDMEC/IST
//    Av. Rovisco Pais, 1, 1049-001 Lisboa, Portugal
//    http://www.dem.ist.utl.pt/idmec/
//
// Based on the initial work of:
//    Lester E. Godwin
//    PushCorp, Inc.
//    Dallas, Texas
//    972-840-0208 x102
//    godwin@pushcorp.com
//    http://www.icsi.berkeley.edu/~storn/devcpp.zip
//
// Last Modified: Oct, 1, 2011
//

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "DESolver.hpp"
#include <errno.h>
#include <ctime>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <mpi.h>

#define EXIT_FAILURE 1

#ifdef _MSC_VER
FILE* file_open_w( const std::string& filename )
{
	FILE* pdump_glb;
	fopen_s( &pdump_glb, filename.c_str(), "w" );
	return pdump_glb;
}
#else
FILE* file_open_w( const std::string& filename )
{
	FILE* pdump_glb = fopen( filename.c_str(), "w" );
	return pdump_glb;
}
#endif

DESolver::DESolver( int dim, int popSize ) :
          nMPI_tasks( 0 ), nMPI_rank( 0 ),
          nDim( dim ), nGlobalPop( popSize ), nLocalPop( 0 ),
          generations( 0 ), DE_Strategy( stRand1Exp ),
          scale( 0.7 ), probability( 0.5 ), bestEnergy( 0.0 ),
          population( popSize, dim ), popEnergy( popSize ),
          trialSolution( dim ), bestSolution( dim ),
          urand( 0 )
{
  time_t rawtime;

#ifndef _MSC_VER
  struct tm* timeinfo;
  time( &rawtime );
  timeinfo = localtime( &rawtime );
  strftime( time_str, 21, "%Y%m%d_%H%M%S", timeinfo );
#else
  struct tm timeinfo;
  time( &rawtime );
  localtime_s( &timeinfo, &rawtime );
  strftime( time_str, 21, "%Y%m%d_%H%M%S", &timeinfo );
#endif
}

void DESolver::MPI_Init( int argc, char* argv[] )
{
  if( ::MPI_Init( &argc, &argv ) != MPI_SUCCESS )
  {
    printf( "MPI_Init failed\n" );
    MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
  }

  MPI_Comm_size( MPI_COMM_WORLD, &nMPI_tasks );
  MPI_Comm_rank( MPI_COMM_WORLD, &nMPI_rank );

//#ifndef ONERUN
//  if( nMPI_tasks < 2 )
//  {
//    fprintf(stderr, "DESolver - parallel code needs at least two tasks\n" );
//    MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
//  }
//#endif

  if( nGlobalPop % nMPI_tasks != 0 )
  {
    fprintf(stderr, "DESolver - population size must be a multiple of the number of tasks\n" );
    MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
  }

  nLocalPop = nGlobalPop / nMPI_tasks;
  urand.seed( ( nMPI_rank + 23 ) * ( nMPI_rank + 19 ) );
}

void DESolver::MPI_Finalize()
{
  ::MPI_Finalize();
}

void DESolver::Setup( Strategy deStrategy, double diffScale, double crossoverProb )
{
  int i;

  DE_Strategy	= deStrategy;
  scale		    = diffScale;
  probability = crossoverProb;

  if( !nMPI_rank )
  {
    for( i = 0; i < nGlobalPop; i++ )
      for( int j = 0; j < nDim; j++ )
        population[i][j] = urand.dbl();

  }

  MPI_Bcast( time_str, 22, MPI_CHAR, 0, MPI_COMM_WORLD );

  MPI_Bcast( population(), nGlobalPop*nDim, MPI_DOUBLE, 0, MPI_COMM_WORLD );

  for( i = 0; i < nGlobalPop; i++ )
    popEnergy[i] = 1.0E20;

  for( i = 0; i < nDim; i++ )
    bestSolution[i] = 0.0;

  switch (DE_Strategy)
  {
    case stBest1Exp:
      calcTrialSolution = &DESolver::Best1Exp;
      break;
    case stRand1Exp:
      calcTrialSolution = &DESolver::Rand1Exp;
      break;
    case stRandToBest1Exp:
      calcTrialSolution = &DESolver::RandToBest1Exp;
      break;
    case stBest2Exp:
      calcTrialSolution = &DESolver::Best2Exp;
      break;
    case stRand2Exp:
      calcTrialSolution = &DESolver::Rand2Exp;
      break;
    case stBest1Bin:
      calcTrialSolution = &DESolver::Best1Bin;
      break;
    case stRand1Bin:
      calcTrialSolution = &DESolver::Rand1Bin;
      break;
    case stRandToBest1Bin:
      calcTrialSolution = &DESolver::RandToBest1Bin;
      break;
    case stBest2Bin:
      calcTrialSolution = &DESolver::Best2Bin;
      break;
    case stRand2Bin:
      calcTrialSolution = &DESolver::Rand2Bin;
      break;
    case stBest3Trig:
      calcTrialSolution = &DESolver::Best3Trig;
      break;
    case stBest1Grad:
      calcTrialSolution = &DESolver::Best1Grad;
      break;
  }
}

void DESolver::limitTrialSolution( Vector& solution, int candidate )
{
  bool found;

  for( int i = 0; i < nDim; i++ )
  {
    found = false;
    do {

      if( solution[i] < 0.0 )
      {
        if( urand.dbl() > 0.618 )
          solution[i] = 0.0;
        else
          solution[i] += 0.5 * urand.dbl();
      }
      else
      if( solution[i] > 1.0 )
      {
        if( urand.dbl() > 0.618 )
          solution[i] = 1.0;
        else
          solution[i] -= 0.5 * urand.dbl();
      }
      else
        found = true;

    } while( !found );
  }
}

bool DESolver::Solve(int maxGenerations, const char* dfile, bool verbose )
{
  int generation;
  int candidate;
  bool bAtSolution;

  int nSend  = nLocalPop * nDim;
  int firstC = nMPI_rank * nLocalPop;
  int lastC  = firstC + nLocalPop;

  bestEnergy = 1.0E20;
  bAtSolution = false;

  double t_start = MPI_Wtime();
  double t_elapsed;

  FILE* pdump_opt = 0;

  if( !nMPI_rank )
  {
    std::string filename_opt( dfile );
    filename_opt += "_";
    filename_opt += time_str;
    filename_opt += ".txt";

    std::cout << "dumpfile optims: " << filename_opt << std::endl << std::flush;

    pdump_opt = file_open_w( filename_opt.c_str() );
    if( pdump_opt == NULL )
    {
      printf( "Dump file open error\n" );
      MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
    }

    fprintf( pdump_opt, "DE strategy: %i\n", DE_Strategy );
    fprintf( pdump_opt, "Dim: %i\n", nDim );
    fprintf( pdump_opt, "maxGenerations: %i\n", maxGenerations );
    fprintf( pdump_opt, "GlobalPopSize: %i\n", nGlobalPop );
    fprintf( pdump_opt, "LocalPopSize: %i\n", nLocalPop );
    fprintf( pdump_opt, "scale: %f\n", scale );
    fprintf( pdump_opt, "probability: %f\n", probability );
    fprintf( pdump_opt, "MPI tasks: %i\n", nMPI_tasks );
    dump_comments( pdump_opt );
    fflush( pdump_opt );
  }

  std::stringstream oss;
  oss << dfile << "_" << time_str << "_proc_" << nMPI_rank << ".txt";
  std::string filename_glb( oss.str() );

  std::cout << "dumpfile global: " << filename_glb.c_str() << std::endl << std::flush;

  FILE* pdump_glb = file_open_w( filename_glb );
  if( pdump_glb == NULL )
  {
    printf( "Dump file open error\n" );
    MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
  }

  for( generation=0; generation < maxGenerations && !bAtSolution; generation++ )
  {
    int newindividuals = 0;

    for( candidate=firstC; candidate < lastC; candidate++ )
    {
      if( generation > 0 )
      {
        (this->*calcTrialSolution)( candidate );
        limitTrialSolution( trialSolution, candidate );
      } else
        trialSolution.copy( population[ candidate ] );

      trialEnergy = EnergyFunction( trialSolution, bAtSolution );

      std::cout << nMPI_rank;
      for( int i = 0; i < nDim; i++ )
      {
        std::cout << std::setprecision(8) << std::setiosflags(std::ios::fixed )
                  << "\t" << trialSolution[i];
        fprintf( pdump_glb, "%f\t", trialSolution[i] );
      }

      std::cout << "\tf(" << std::setw(3) << candidate << ") = "
                << trialEnergy << std::endl << std::flush;
      fprintf( pdump_glb, "\t%f\n", trialEnergy );
      fflush( pdump_glb );

      if( trialEnergy < popEnergy[candidate] )
      {
        popEnergy[candidate] = trialEnergy;
        population.rowCopy( candidate, trialSolution() );

        if( trialEnergy < bestEnergy)
        {
          bestEnergy = trialEnergy;
          bestSolution.copy( trialSolution() );
        }

        newindividuals++;
      }

      MPI_Allgather( population() + nMPI_rank * nSend, nSend, MPI_DOUBLE,
                     population(), nSend, MPI_DOUBLE,
                     MPI_COMM_WORLD );

      MPI_Allgather( popEnergy() + nMPI_rank * nLocalPop, nLocalPop, MPI_DOUBLE,
                     popEnergy(), nLocalPop, MPI_DOUBLE,
                     MPI_COMM_WORLD );

      for( int c=0; c< nGlobalPop; c++ )
      {
        if( popEnergy[c] < bestEnergy)
        {
          bestEnergy = popEnergy[c];
          bestSolution.copy( population[c] );
        }
      }
    }

    t_elapsed = MPI_Wtime() - t_start;

    if( !nMPI_rank && verbose)
    {
      std::cout << "\nGeneration: " << generation << std::endl;
      std::cout << "Best Coefficients:"  << std::endl;
      dump_solution( stdout, bestSolution );
      std::cout << "Function value: " << std::setprecision(9) << bestEnergy  << std::endl;
      std::cout << "New individuals: " << newindividuals << std::endl;
      std::cout << "Time elapsed: " << std::setprecision(1) << t_elapsed << std::endl << std::flush;
    }

    if( !nMPI_rank )
    {
      fprintf( pdump_opt, "\nGeneration: %i\n", generation );
      fprintf( pdump_opt, "Best Coefficients:\n");
      dump_solution( pdump_opt, bestSolution );
      fprintf( pdump_opt, "Function value = %.9e\n", bestEnergy );
      fprintf( pdump_opt, "New individuals = %i\n", newindividuals );
      fprintf( pdump_opt, "Time elapsed = %f\n", t_elapsed );
      fflush( pdump_opt );
    }
  }

  if( !nMPI_rank )
  {
    fclose( pdump_opt );
  }
  fclose( pdump_glb );

  generations = generation;

  return bAtSolution;
}

void DESolver::Best1Exp( int candidate )
{
  int r1, r2;

  SelectSamples( candidate, &r1, &r2 );
  int n = int( urand.dbl( 0.0, double(nDim) ) );

  trialSolution.copy( population[candidate] );
  for( int i = 0; urand.dbl() < probability && i < nDim; i++)
  {
    if( urand.dbl() < probability )
      trialSolution[n] = bestSolution[n]
                + scale * ( population[r1][n] - population[r2][n] );
    n = (n + 1) % nDim;
  }
}

void DESolver::Best1Grad( int candidate )
{
  int r1;

  SelectSamples( candidate, &r1 );
  int n = int( urand.dbl( 0.0, double(nDim) ) );

  trialSolution.copy( population[candidate] );
  for( int i = 0; urand.dbl() < probability && i < nDim; i++)
  {
    if( urand.dbl( 0.0, 1.0 ) < probability )
      trialSolution[n] = bestSolution[n]
                + scale * ( bestSolution[n] - population[r1][n] );
    n = (n + 1) % nDim;
  }
}

void DESolver::Rand1Exp( int candidate )
{
  int r1, r2, r3;

  SelectSamples( candidate, &r1, &r2, &r3 );
  int n = int( urand.dbl( 0.0, double(nDim) ) );

  trialSolution.copy( population[candidate] );
  for( int i = 0; urand.dbl() < probability && i < nDim; i++)
  {
    if( urand.dbl( 0.0, 1.0 ) < probability )
      trialSolution[n] = population[r1][n]
                       + scale * ( population[r2][n] - population[r3][n] );
    n = (n + 1) % nDim;
  }
}

void DESolver::RandToBest1Exp( int candidate )
{
  int r1, r2;

  SelectSamples( candidate, &r1, &r2 );
  int n = int( urand.dbl( 0.0, double(nDim) ) );

  trialSolution.copy( population[candidate] );
  for( int i = 0; urand.dbl() < probability && i < nDim; i++)
  {
    if(urand.dbl( 0.0, 1.0 ) < probability)
      trialSolution[n] += scale * ( bestSolution[n] - trialSolution[n] )
                        + scale * ( population[r1][n] - population[r2][n] );
    n = (n + 1) % nDim;
  }
}

void DESolver::Best2Exp( int candidate )
{
  int r1, r2, r3, r4;

  SelectSamples( candidate, &r1, &r2, &r3, &r4 );
  int n = int( urand.dbl( 0.0, double(nDim) ) );

  trialSolution.copy( population[candidate] );
  for( int i = 0; urand.dbl() < probability && i < nDim; i++)
  {
    if( urand.dbl( 0.0, 1.0 ) < probability )
      trialSolution[n] = bestSolution[n] +
                    scale * ( population[r1][n] + population[r2][n]
                            - population[r3][n] - population[r4][n] );
    n = (n + 1) % nDim;
  }
}

void DESolver::Rand2Exp( int candidate )
{
  int r1, r2, r3, r4, r5;

  SelectSamples( candidate, &r1, &r2, &r3, &r4, &r5 );
  int n = int( urand.dbl( 0.0, double(nDim) ) );

  trialSolution.copy( population[candidate] );
  for( int i = 0; urand.dbl() < probability && i < nDim; i++)
  {
    if( urand.dbl( 0.0, 1.0 ) < probability )
      trialSolution[n] = population[r1][n]
                       + scale * ( population[r2][n] + population[r3][n]
                                 - population[r4][n] - population[r5][n] );
    n = (n + 1) % nDim;
  }
}

void DESolver::Best1Bin( int candidate )
{
  int r1, r2;

  SelectSamples( candidate, &r1, &r2 );
  int n = int( urand.dbl( 0.0, double(nDim) ) );

  trialSolution.copy( population[candidate] );
  for( int i = 0; urand.dbl() < probability && i < nDim; i++)
  {
    if( urand.dbl( 0.0, 1.0 ) < probability || i == (nDim - 1) )
      trialSolution[n] = bestSolution[n]
                       + scale * ( population[r1][n] - population[r2][n] );
    n = (n + 1) % nDim;
  }
}

void DESolver::Rand1Bin( int candidate )
{
  int r1, r2, r3;

  SelectSamples( candidate, &r1, &r2, &r3 );
  int n = int( urand.dbl( 0.0, double(nDim) ) );

  trialSolution.copy( population[candidate] );
  for( int i = 0; urand.dbl() < probability && i < nDim; i++)
  {
    if( urand.dbl( 0.0, 1.0 ) < probability || i == (nDim - 1) )
      trialSolution[n] = population[r1][n]
                       + scale * ( population[r2][n] - population[r3][n] );
    n = (n + 1) % nDim;
  }
}

void DESolver::RandToBest1Bin( int candidate )
{
  int r1, r2;

  SelectSamples( candidate, &r1, &r2 );
  int n = int( urand.dbl( 0.0, double(nDim) ) );

  trialSolution.copy( population[candidate] );
  for( int i = 0; urand.dbl() < probability && i < nDim; i++)
  {
    if( urand.dbl( 0.0, 1.0 ) < probability || i == (nDim - 1) )
      trialSolution[n] += scale * ( bestSolution[n] - trialSolution[n] )
                        + scale * ( population[r1][n] - population[r2][n] );
    n = (n + 1) % nDim;
  }
}

void DESolver::Best2Bin( int candidate )
{
  int r1, r2, r3, r4;

  SelectSamples( candidate, &r1, &r2, &r3, &r4 );
  int n = int( urand.dbl( 0.0, double(nDim) ) );

  trialSolution.copy( population[candidate] );
  for( int i = 0; urand.dbl() < probability && i < nDim; i++)
  {
    if( urand.dbl( 0.0, 1.0 ) < probability || i == (nDim - 1) )
      trialSolution[n] = bestSolution[n]
                + scale * ( population[r1][n] + population[r2][n]
                          - population[r3][n] - population[r4][n] );
    n = (n + 1) % nDim;
  }
}

void DESolver::Rand2Bin( int candidate )
{
  int r1, r2, r3, r4, r5;

  SelectSamples( candidate, &r1, &r2, &r3, &r4, &r5 );
  int n = int( urand.dbl( 0.0, double(nDim) ) );

  trialSolution.copy( population[candidate] );
  for( int i = 0; urand.dbl() < probability && i < nDim; i++)
  {
    if( urand.dbl( 0.0, 1.0 ) < probability || i == (nDim - 1) )
      trialSolution[n] = population[r1][n]
                       + scale * ( population[r2][n] + population[r3][n]
                                 - population[r4][n] - population[r5][n] );
    n = (n + 1) % nDim;
  }
}

void DESolver::Best3Trig( int candidate )
{
  int r1, r2, r3;

  SelectSamples( candidate, &r1, &r2, &r3 );
  int n = int( urand.dbl( 0.0, double(nDim) ) );

  double p1 = fabs( popEnergy[r1] );
  double p2 = fabs( popEnergy[r2] );
  double p3 = fabs( popEnergy[r3] );
  double pl = p1 + p2 + p3;

  trialSolution.copy( population[candidate] );
  for( int i = 0; urand.dbl() < probability && i < nDim; i++)
  {
    if( urand.dbl( 0.0, 1.0 ) < probability )
      trialSolution[n] =
          ( population[r1][n] + population[r2][n] + population[r3][n] ) / 3.0
             + (p2 - p1) / pl * ( population[r2][n] + population[r1][n] )
             + (p3 - p2) / pl * ( population[r3][n] + population[r2][n] )
             + (p1 - p3) / pl * ( population[r1][n] + population[r3][n] );
    n = (n + 1) % nDim;
  }
}

void DESolver::SelectSamples( int candidate, int *r1, int *r2,
                    int *r3, int *r4, int *r5 )
{
  if( r1 )
    do *r1 = int( urand.dbl( 0.0, double(nGlobalPop) ) );
    while( *r1 == candidate );

  if( r2 )
    do *r2 = int( urand.dbl( 0.0, double(nGlobalPop) ) );
    while( *r2 == candidate || *r2 == *r1 );

  if( r3 )
    do *r3 = int( urand.dbl( 0.0, double(nGlobalPop) ) );
    while( *r3 == candidate || *r3 == *r2 || *r3 == *r1 );

  if( r4 )
    do *r4 = int( urand.dbl( 0.0, double(nGlobalPop) ) );
    while( *r4 == candidate || *r4 == *r3 || *r4 == *r2 || *r4 == *r1 );

  if( r5 )
    do *r5 = int( urand.dbl( 0.0, double(nGlobalPop) ) );
    while( *r5 == candidate || *r5 == *r4 || *r5 == *r3 || *r5 == *r2 || *r5 == *r1 );
}
