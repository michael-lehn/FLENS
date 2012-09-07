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
    using std::abs;

    SDenseVectorView  x(SArrayView(*N, X, abs(*INCX)), *INCX<0);
    SDenseVectorView  y(SArrayView(*N, Y, abs(*INCY)), *INCY<0);

    blas::swap(x, y);
}


void
BLAS(dswap)(const INTEGER   *N,
            double          *X,
            const INTEGER   *INCX,
            double          *Y,
            const INTEGER   *INCY)
{
    using std::abs;

    DDenseVectorView  x(DArrayView(*N, X, abs(*INCX)), *INCX<0);
    DDenseVectorView  y(DArrayView(*N, Y, abs(*INCY)), *INCY<0);

    blas::swap(x, y);
}

void
BLAS(cswap)(const INTEGER   *N,
            cfloat          *X,
            const INTEGER   *INCX,
            cfloat          *Y,
            const INTEGER   *INCY)
{
    using std::abs;

    CDenseVectorView  x(CArrayView(*N, X, abs(*INCX)), *INCX<0);
    CDenseVectorView  y(CArrayView(*N, Y, abs(*INCY)), *INCY<0);

    blas::swap(x, y);
}


void
BLAS(zswap)(const INTEGER   *N,
            cdouble         *X,
            const INTEGER   *INCX,
            cdouble         *Y,
            const INTEGER   *INCY)
{
    using std::abs;

    ZDenseVectorView  x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);
    ZDenseVectorView  y(ZArrayView(*N, Y, abs(*INCY)), *INCY<0);

    blas::swap(x, y);
}

} // extern "C"
