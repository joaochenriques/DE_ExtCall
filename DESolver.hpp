//////////////////////////////////////////////////////////////////////
//
// Differential Evolution Test Program
// Based on algorithms developed by Dr. Rainer Storn & Kenneth Price
//
// Parallel version By:
//    Joao C. C. Henriques
//    joaochenriques@ist.utl.pt, joaochenriques@gmail.com
//    Institute of Mechanical Engineering / Instituto Superior Tecnico
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
#if !defined(_DESOLVER_HPP)
#define _DESOLVER_H

#include "nr_urand.hpp"
#include "matvec.hpp"
#include <stdio.h>
#include <string>

class DESolver
{
public:
  typedef enum {
    stBest1Exp = 1, stRand1Exp, stRandToBest1Exp, stBest2Exp, stRand2Exp,
    stBest1Bin, stRand1Bin, stRandToBest1Bin, stBest2Bin, stRand2Bin,
    stBest3Trig, stBest1Grad
  } Strategy;

  DESolver( int dim, int popSize );
  virtual ~DESolver() {}

  // Setup() must be called before solve to set min, max, strategy etc.
  virtual void Setup( Strategy deStrategy, double diffScale, double crossoverProb );

  // Solve() returns true if EnergyFunction() returns true.
  // Otherwise it runs maxGenerations generations and returns false.
  bool Solve( int maxGenerations, const char* dfile, bool verbose );

  // EnergyFunction must be overridden for problem to solve
  // testSolution[] is nDim array for a candidate solution
  // setting bAtSolution = true indicates solution is found
  // and Solve() immediately returns true.
  virtual double EnergyFunction( Vector& testSolution,
                                 bool &bAtSolution ) { return 1E20; }

  virtual void dump_comments( FILE* ) {}
  virtual void dump_solution( FILE*, const Vector& ) {}

  int Dimension() { return nDim; }
  int Population() { return nGlobalPop; }

  // Call these functions after Solve() to get results.
  double Energy() { return bestEnergy; }
  const Vector& Solution() { return bestSolution; }

  int Generations() { return generations; }

  int MPI_Rank() { return nMPI_rank; }
  int MPI_Tasks() { return nMPI_tasks; }

  void MPI_Init( int argc, char* argv[] );
  void MPI_Finalize();

  void dump_energy( FILE* );

protected:
  typedef void ( DESolver::*StrategyFunction )( int );

  void SelectSamples( int candidate, int *r1, int *r2 = 0,
                      int *r3 = 0, int *r4 = 0, int *r5 = 0 );

  void limitTrialSolution( Vector& solution, int candidate );

  int nMPI_tasks;
  int nMPI_rank;

  int nDim;
  int nGlobalPop;
  int nLocalPop;
  int generations;

  Strategy DE_Strategy;
  StrategyFunction calcTrialSolution;

  double scale;
  double probability;

  double bestEnergy;
  double trialEnergy;

  Matrix population;
  Vector popEnergy;
  Vector trialSolution;
  Vector bestSolution;

  URand urand;
  char time_str[22];

private:
  void Best1Grad( int candidate );
  void Best1Exp( int candidate );
  void Rand1Exp( int candidate );
  void RandToBest1Exp( int candidate );
  void Best2Exp( int candidate );
  void Rand2Exp( int candidate );
  void Best1Bin( int candidate );
  void Rand1Bin( int candidate );
  void RandToBest1Bin( int candidate );
  void Best2Bin( int candidate );
  void Rand2Bin( int candidate );
  void Best3Trig( int candidate );
};

#endif // _DESOLVER_HPP
