#include <assert.h>

class URand {
protected:
  unsigned long long  u,v,w;

public:
  URand(unsigned long long j)
  {
    assert( sizeof(long long) >= 8 );
    seed(j);
  }

  void seed(unsigned long long j)
  {
    v = 4101842887655102017LL;
    w = 1;
    u = j^v; int64();
    v = u; int64();
    w = v; int64();
  }

  inline unsigned long long int64()
  {
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    w = 4294957665U*(w & 0xffffffff) + (w >> 32);
    unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
    return (x + v) ^ w;
  }

  // Return random double-precision floating value in the range 0. to 1.
  inline double dbl()
  {
    return 5.42101086242752217E-20 * int64();
  }

  inline double dbl( double min, double max )
  {
    return min + dbl() * ( max - min );
  }

  //inline unsigned int int32()
  //{
  //  return unsigned int( int64() );
  //}
};

