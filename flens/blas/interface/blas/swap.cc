#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(sswap)(const INTEGER   *N,
            float           *X,
            const INTEGER   *INCX,
            float           *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: sswap");

#   ifdef TEST_DIRECT_CBLAS

    cblas_sswap(*N, X, *INCX, Y, *INCY);

#   else

        using std::abs;

        SDenseVectorView  x(SArrayView(*N, X, abs(*INCX)), *INCX<0);
        SDenseVectorView  y(SArrayView(*N, Y, abs(*INCY)), *INCY<0);

        blas::swap(x, y);
#   endif
}


void
BLAS(dswap)(const INTEGER   *N,
            double          *X,
            const INTEGER   *INCX,
            double          *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: dswap");

#   ifdef TEST_DIRECT_CBLAS

        cblas_dswap(*N, X, *INCX, Y, *INCY);

#   else

        using std::abs;

        DDenseVectorView  x(DArrayView(*N, X, abs(*INCX)), *INCX<0);
        DDenseVectorView  y(DArrayView(*N, Y, abs(*INCY)), *INCY<0);

        blas::swap(x, y);
#   endif
}

void
BLAS(cswap)(const INTEGER   *N,
            cfloat          *X,
            const INTEGER   *INCX,
            cfloat          *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: cswap");

#   ifdef TEST_DIRECT_CBLAS

        cblas_cswap(*N,
                    reinterpret_cast<float *>(X), *INCX,
                    reinterpret_cast<float *>(Y), *INCY);

#   else

        using std::abs;

        CDenseVectorView  x(CArrayView(*N, X, abs(*INCX)), *INCX<0);
        CDenseVectorView  y(CArrayView(*N, Y, abs(*INCY)), *INCY<0);

        blas::swap(x, y);
#   endif
}


void
BLAS(zswap)(const INTEGER   *N,
            cdouble         *X,
            const INTEGER   *INCX,
            cdouble         *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: zswap");

#   ifdef TEST_DIRECT_CBLAS

        cblas_zswap(*N,
                    reinterpret_cast<double *>(X), *INCX,
                    reinterpret_cast<double *>(Y), *INCY);

#   else

        using std::abs;

        ZDenseVectorView  x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);
        ZDenseVectorView  y(ZArrayView(*N, Y, abs(*INCY)), *INCY<0);

        blas::swap(x, y);
#   endif
}

} // extern "C"
