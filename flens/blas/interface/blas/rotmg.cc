#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(srotmg)(float   *A,
             float   *B,
             float   *C,
             float  D,
             float   *S)
{
    ASSERT(0);
}

void
BLAS(drotmg)(double  *A,
             double  *B,
             double  *C,
             double  D,
             double  *S)
{
    ASSERT(0);
}

} // extern "C"
