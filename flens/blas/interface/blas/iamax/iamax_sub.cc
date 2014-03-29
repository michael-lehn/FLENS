#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(isamax_sub)(const INTEGER   *N,
                 const float     *X,
                 const INTEGER   *INCX,
                 INTEGER         *RES)
{
#   ifdef TEST_DIRECT_CBLAS

        *RES = cblas_isamax(*N, X, *INCX);

#   else

        using std::abs;

        SDenseVectorConstView  x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);

        *RES = blas::iamax(x);

#   endif
}

void
BLAS(idamax_sub)(const INTEGER   *N,
                 const double    *X,
                 const INTEGER   *INCX,
                 INTEGER         *RES)
{
#   ifdef TEST_DIRECT_CBLAS

        *RES = cblas_idamax(*N, X, *INCX);

#   else

        using std::abs;

        DDenseVectorConstView  x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);

        *RES = blas::iamax(x);

#   endif
}

void
BLAS(icamax_sub)(const INTEGER   *N,
                 const cfloat    *X,
                 const INTEGER   *INCX,
                 INTEGER         *RES)
{
#   ifdef TEST_DIRECT_CBLAS

        *RES = cblas_icamax(*N, reinterpret_cast<const float *>(X), *INCX);

#   else

        using std::abs;

        CDenseVectorConstView  x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);

        *RES = blas::iamax(x);
#   endif
}

void
BLAS(izamax_sub)(const INTEGER   *N,
                 const cdouble   *X,
                 const INTEGER   *INCX,
                 INTEGER         *RES)
{
#   ifdef TEST_DIRECT_CBLAS

        *RES = cblas_izamax(*N, reinterpret_cast<const double *>(X), *INCX);

#   else

        using std::abs;

        ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);

        *RES = blas::iamax(x);

#   endif
}

} // extern "C"
