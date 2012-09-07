#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(srotg)(float   *A,
            float   *B,
            float   *C,
            float   *S)
{
    blas::rotg(*A, *B, *C, *S);
}

void
BLAS(drotg)(double  *A,
            double  *B,
            double  *C,
            double  *S)
{
    blas::rotg(*A, *B, *C, *S);
}

void
BLAS(crotg)(cfloat   *A,
            cfloat   *B,
            float    *C,
            cfloat   *S)
{
    blas::rotg(*A, *B, *C, *S);
}

void
BLAS(zrotg)(cdouble  *A,
            cdouble  *B,
            double   *C,
            cdouble  *S)
{
    blas::rotg(*A, *B, *C, *S);
}


} // extern "C"
