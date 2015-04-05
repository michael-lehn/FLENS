#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(srotg)(float   *A,
            float   *B,
            float   *C,
            float   *S)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: srotg");

#   ifdef TEST_DIRECT_CBLAS

        cblas_srotg(A, B, C, S);

#   else

        blas::rotg(*A, *B, *C, *S);

#   endif
}

void
BLAS(drotg)(double  *A,
            double  *B,
            double  *C,
            double  *S)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: drotg");

#   ifdef TEST_DIRECT_CBLAS

        cblas_drotg(A, B, C, S);

#   else

        blas::rotg(*A, *B, *C, *S);

#   endif
}

void
BLAS(crotg)(cfloat   *A,
            cfloat   *B,
            float    *C,
            cfloat   *S)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: crotg");

    blas::rotg(*A, *B, *C, *S);
}

void
BLAS(zrotg)(cdouble  *A,
            cdouble  *B,
            double   *C,
            cdouble  *S)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: drotg");

    blas::rotg(*A, *B, *C, *S);
}


} // extern "C"
