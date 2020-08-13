#if !defined(_MATVEC_HPP)
#define _MATVEC_HPP

#include <assert.h>

class Vector 
{
protected:
  double* buffer;
  int dim;

public:
  int size() const { return dim; }
  
  void resize( int d )
  {
    if( d != dim )
    {
      dim = d;
      delete[] buffer;
      buffer = new double[ dim ];
    }
  }

  Vector( int d ) : dim(d)
  {
    buffer = new double[ dim ];
  }
    
  Vector( int d, double c ) : dim(d)
  {
    buffer = new double[ dim ];
    for( int j = 0; j < dim; j++ )
      buffer[j] = c;
  }
    
  Vector( int d, const double* c ) : dim(d)
  {
    buffer = new double[ dim ];
    for( int j = 0; j < dim; j++ )
      buffer[j] = c[j];
  }
    
  Vector( const Vector& v ) : dim( v.dim )
  {
    buffer = new double[ dim ];
    copy( v.buffer ); 
  }
  
  ~Vector()
  {
    delete buffer;
  }
  
  inline double& operator[]( int i )
  {
    return buffer[i];
  }
  
  inline double operator()( int i ) const
  {
    return buffer[i];
  }
  
  inline double* operator() () { return buffer; }
  
  void copy( double* b )
  {
    for( int j = 0; j < dim; j++ )
      buffer[j] = b[j];
    //memcpy( buffer, b, dim * sizeof(double) );
  }

  void operator = ( const Vector& b )
  {
    if( b.dim != dim )
      resize( b.dim );
    for( int j = 0; j < dim; j++ )
      buffer[j] = b.buffer[j];
  }
};

class Matrix 
{
protected:
  double* buffer;
  int rows;
  int cols;

public:

  void size( int& r, int& c )
  {
    r = rows;
    c = cols;
  }
  
  Matrix( int r, int c ) : rows(r), cols(c)
  {
    buffer = new double[ rows * cols ];
  }
    
  ~Matrix()
  {
    delete buffer;
  }
  
  inline double* operator[]( int i )
  {
    return buffer + i*cols;
  }
  
  inline double* operator() () { return buffer; }

  void rowCopy( int i, double* b )
  {
    double *a = buffer + i*cols;
    for( int j = 0; j < cols; j++ ) a[j] = b[j];
    //memcpy( a, b, cols * sizeof(double) );
  }
};

#endif

