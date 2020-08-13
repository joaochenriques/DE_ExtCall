////////////////////////////////////////////////////////////////////////////////
//
// Differential Evolution Test Program
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
//  mpicxx DE_ExtCall.cpp DESolver.cpp -lconfig++ -o DE-opt
//
//******************************************************************************
// Linux remote cmd line with 6 cpu: nohup mpirun -n 6 DE-opt &

//  nohup mpirun -n 6 DE-opt nohup.out & tail -f nohup.out
//******************************************************************************
// LIB =
//  C:\Program Files\Microsoft Visual Studio 8\VC\lib\amd64;
//  C:\Program Files\Microsoft Visual Studio 8\VC\lib;
//  C:\Program Files\Microsoft.NET\SDK\v2.0 64bit\Lib;
//
// INCLUDE =
//  C:\Program Files (x86)\Microsoft Visual Studio 8\VC\include;
//  C:\Program Files\Microsoft.NET\SDK\v2.0 64bit\include;
//  C:\Program Files (x86)\Microsoft Visual Studio 8\VC\PlatformSDK\Include
//
// PATH +=
//  C:\Program Files\Microsoft HPC Pack 2008 SDK\Bin\;
//  C:\Program Files\Microsoft.NET\SDK\v2.0 64bit\Bin;
//  C:\Windows\Microsoft.NET\Framework64\v2.0.50727;
//  C:\Program Files\Microsoft Visual Studio 8\VC\bin;
//  C:\Program Files\Microsoft Visual Studio 8\Common7\IDE;
//  C:\Program Files\Microsoft Visual Studio 8\VC\vcpackages;
//******************************************************************************
//
#include <iostream>
#include <sstream>
#include <cstdio>
#include <unistd.h>
#include <iomanip>
#include <fstream>
#include <vector>
#include <iterator>

#include "DESolver.hpp"
#include <libconfig.h++>

using namespace libconfig;

#ifdef _MSC_VER
#include <windows.h>
std::wstring s2ws(const std::string& s)
{
    int len;
    int slength = (int)s.length() + 1;
    len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0);
    wchar_t* buf = new wchar_t[len];
    MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
    std::wstring r(buf);
    delete[] buf;
    return r;
}

int call_exec( const std::string& cmd, const std::string& args )
{
	std::wstring wcmd = s2ws( cmd );
	std::wstring wargs = s2ws( args );

	STARTUPINFO StartInfo;
	PROCESS_INFORMATION ProcInfo;

	memset(&ProcInfo, 0, sizeof(ProcInfo));
	memset(&StartInfo, 0 , sizeof(StartInfo) );
	StartInfo.cb = sizeof(StartInfo);
    StartInfo.dwFlags = STARTF_USESHOWWINDOW;
	int res = CreateProcess( wcmd.c_str(), const_cast<LPWSTR>( wargs.c_str() ),
                           NULL, NULL, FALSE,
                           CREATE_BREAKAWAY_FROM_JOB | NORMAL_PRIORITY_CLASS,
                           NULL, NULL, &StartInfo, &ProcInfo);

	if (res)
    {
		WaitForSingleObject(ProcInfo.hThread, INFINITE);
	  return 0;
    } else
    return 1;
}
#else
#include <cstdlib>
int call_exec( const std::string& cmd, const std::string& args )
{
  std::string cmd_line( cmd );
  cmd_line += args;
  return std::system( cmd_line.c_str() );
}
#endif

//~############################################################################
class Params {
protected:
  Config cfg;

public:

  Params( const char* filename )
  {
    // Read the file. If there is an error, report it and exit.
    try
    {
      cfg.readFile( filename );
    }

    catch(const FileIOException &fioex)
    {
      std::cerr << "I/O error while reading file." << std::endl;
      exit(1);
    }
    catch(const ParseException &pex)
    {
      std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                << " - " << pex.getError() << std::endl;
      exit(1);
    }
  }

  double get_double( const std::string& section_name, const std::string& name )
  {
    const Setting& root = cfg.getRoot();
    double val;
    try
    {
      Setting &section = root[ section_name.c_str() ];
      if( !section.lookupValue( name.c_str(), val ) )
      {
        std::cout << "ERROR@Params: name \"" << name
                  << "\" not found in section \""
                  << section_name << "\"\n";
        exit(1);
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      std::cout << "ERROR@Params: SettingNotFoundException in section \""
           << section_name << "\"\n";
      exit(1);
    }
    return val;
  }

  int get_int( const std::string& section_name, const std::string& name )
  {
    const Setting& root = cfg.getRoot();
    int val;
    try
    {
      const Setting &section = root[ section_name.c_str() ];
      if( !section.lookupValue( name.c_str(), val ) )
      {
        std::cout << "ERROR@Params: name \"" << name
                  << "\" not found in section \""
                  << section_name << "\"\n";
        exit(1);
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      std::cout << "ERROR@Params: SettingNotFoundException in section \""
           << section_name << "\"\n";
      exit(1);
    }
    return val;
  }

  std::string get_string( const std::string& section_name,
                          const std::string& name )
  {
    const Setting& root = cfg.getRoot();
    std::string val;
    try
    {
      const Setting &section = root[ section_name.c_str() ];
      if( !section.lookupValue( name.c_str(), val ) )
      {
        std::cout << "ERROR@Params: name \"" << name
                  << "\" not found in section \""
                  << section_name << "\"\n";
        exit(1);
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      std::cout << "ERROR@Params: SettingNotFoundException in section \""
           << section_name << "\"\n";
      exit(1);
    }
    return val;
  }
};

//~############################################################################
void ExitError( const std::string& msg )
{
  std::cerr << msg << std::endl;
  exit(1);
}

//~############################################################################
bool get_tokens( std::stringstream& buffer, std::vector<std::string>& str_lst )
{
  std::string line;

  if( getline( buffer, line ) )
  {
    std::stringstream iss;
    iss << line;

    std::istream_iterator<std::string> line_begin( iss );
    std::istream_iterator<std::string> line_end;

    str_lst.assign( line_begin, line_end );
    return false;
  }
  return true;
}

//~############################################################################
template <typename T>
bool convert_string( const std::string& s, T& value )
{
  std::istringstream is(s);
  if( !(is >> value) || !is.eof() ) return true;
  return false;
}

//~############################################################################
class DE_Biradial : public DESolver
{
protected:
  Vector Xmin;
  Vector Xmax;
  Vector DeltaX;

  double InvalidValue;

  std::string RunName;
  std::string ExtCall;
  std::string ExtArgs;

  std::string file_name;
  std::string file_name_cout;
  std::string file_name_cerr;
  std::string file_name_cfg;
  std::string file_name_val;

public:

  DE_Biradial(int dim,int pop, Params& params ) :
      DESolver(dim,pop), Xmin(dim), Xmax(dim), DeltaX(dim)
  {
    InvalidValue = 1E6;

    if( nDim <= 0 or nDim >= 100 )
    {
      std::cout << "DE_ExtCall implies: 0 < nDim < 100" << std::endl;
      exit(1);
    }

    ExtCall = params.get_string( "USER", "ExtCall" );
    RunName = params.get_string( "USER", "RunName" );
    InvalidValue = params.get_double( "USER", "InvalidValue" );

    for( int i = 0; i < nDim; i++ )
    {
      std::stringstream xmin_name;
      xmin_name << "VAR_X" << std::setfill('0') << std::setw(2) << i << "_MIN";

      std::stringstream xmax_name;
      xmax_name << "VAR_X" << std::setfill('0') << std::setw(2) << i << "_MAX";

      Xmin[i] = params.get_double( "USER", xmin_name.str() );
      Xmax[i] = params.get_double( "USER", xmax_name.str() );
    }

    for( int i = 0; i < nDim; i++ )
      DeltaX[i] = Xmax[i] - Xmin[i];
  }

  void Setup( Strategy deStrategy, double diffScale, double crossoverProb )
  {
    DESolver::Setup( deStrategy, diffScale, crossoverProb );

    std::stringstream name;
    name << RunName << "_" << time_str << "_"
         << std::setfill('0') << std::setw(3) << nMPI_rank;
    file_name = name.str();

    file_name_cout  = file_name;
    file_name_cout += ".cout";

    file_name_cerr  = file_name;
    file_name_cerr += ".cerr";

    file_name_cfg  = file_name;
    file_name_cfg += ".cfg";

    file_name_val  = file_name;
    file_name_val += ".value";

    ExtArgs  = " ";
    ExtArgs += file_name_cfg;
    ExtArgs += " 1> ";
    ExtArgs += file_name_cout;
    ExtArgs += " 2> ";
    ExtArgs += file_name_cerr;
  }

  double EnergyFunction( Vector& trial,bool &bAtSolution);
  virtual void dump_comments( FILE* pdump );
  virtual void dump_solution( FILE* pdump, const Vector& trial );
  double read_val_file( const std::string filename );
};

void DE_Biradial::dump_solution( FILE* pdump, const Vector& trial )
{
  double f;

  fprintf( pdump, "DE space\n" );
  for( int i = 0; i < trial.size(); i++ )
    fprintf( pdump, "X[%d] = %lf\n", i, trial(i) );

  fprintf( pdump, "Function space\n" );
  for( int i = 0; i < trial.size(); i++ )
  {
    f = Xmin[i] + DeltaX[i] * trial(i);
    fprintf( pdump, "x[%d] = %lf\n", i, f );
  }
}

void DE_Biradial::dump_comments( FILE* pdump )
{
  fprintf( pdump, "DE_ExtCall\n\n" );

  for( int i = 0; i < Xmin.size(); i++ )
  {
    fprintf( pdump, "Xmin[%i] = %f\n", i, Xmin[i] );
    fprintf( pdump, "Xmax[%i] = %f\n", i, Xmax[i] );
  }

  fprintf( pdump, "\n" );
}

double DE_Biradial::EnergyFunction( Vector& trial,bool &bAtSolution)
{
  Vector X( trial.size() );
  for( int i = 0; i < trial.size(); i++ )
    X[i] = Xmin[i] + DeltaX[i] * trial(i);

  FILE* pfile = fopen( file_name_cfg.c_str(), "w" );
  fprintf( pfile, "%i\n", nDim );
  for( int var = 0; var < nDim; var ++ )
    fprintf( pfile, "%.15E\n", X[var] );
  fprintf( pfile, "%s\n", file_name_val.c_str() );
  fclose( pfile );

  int res1 = call_exec( ExtCall, ExtArgs );

  if( res1 != 0 )
    return InvalidValue;
  else
    return read_val_file( file_name_val );
}

double DE_Biradial::read_val_file( const std::string filename )
{
  std::stringstream buffer;
  std::ifstream t( filename.c_str() );
  buffer << t.rdbuf();
  t.close();

  if( buffer )
  {
    double value;
    std::vector<std::string> tokens;
    get_tokens( buffer, tokens );

    if( tokens.size() != 1 )
      ExitError( "ERROR@read_val_file: invalid file data" );

    if( convert_string( tokens[0], value ) )
      ExitError( "ERROR@read_val_file: invalid double" );
    else
      return value;
  }
  return InvalidValue;
}

int main( int argc, char *argv[] )
{
  Params params( "./DE_params.cfg" );

  int strat = params.get_int( "DE", "strategy" );
  int N_DIM = params.get_int( "DE", "N_DIM" );
  int N_POP = params.get_int( "DE", "N_POP" );

  int MAX_GENERATIONS	= params.get_int( "DE", "MAX_GENERATIONS" );

  double diffScale = params.get_double( "DE", "diffScale" );
  double crossoverProb = params.get_double( "DE", "crossoverProb" );

  DE_Biradial solver(N_DIM,N_POP,params);
  solver.MPI_Init( argc, argv );

  if( strat < DESolver::stBest1Exp or strat > DESolver::stBest1Grad )
    ExitError( "ERROR@invalid strategy." );

  solver.Setup( (DESolver::Strategy)strat, diffScale, crossoverProb );
  solver.Solve( MAX_GENERATIONS, "./dumps/dumpfile", true );

  if( !solver.MPI_Rank() )
  {
    bool dummy;
    Vector solution( solver.Solution() );
    printf("\nGeneration: %i\n", solver.Generations() );
    printf("Best Coefficients:\n");
    solver.dump_solution( stdout, solution );
    printf( "Function = %.9e\n\n", solver.EnergyFunction(solution,dummy) );
  }

  solver.MPI_Finalize();
  return 0;
}
