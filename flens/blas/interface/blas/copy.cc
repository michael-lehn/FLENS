#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(scopy)(const INTEGER   *N,
            const float     *X,
            const INTEGER   *INCX,
            float           *Y,
            const INTEGER   *INCY)
{
    using std::abs;

    SDenseVectorConstView  x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    SDenseVectorView       y(SArrayView(*N, Y, abs(*INCY)), *INCY<0);

    y = x;
}


void
BLAS(dcopy)(const INTEGER   *N,
            const double    *X,
            const INTEGER   *INCX,
            double          *Y,
            const INTEGER   *INCY)
{
    using std::abs;

    DDenseVectorConstView  x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    DDenseVectorView       y(DArrayView(*N, Y, abs(*INCY)), *INCY<0);

    y = x;
}

void
BLAS(ccopy)(const INTEGER   *N,
            const cfloat    *X,
            const INTEGER   *INCX,
            cfloat          *Y,
            const INTEGER   *INCY)
{
    using std::abs;

    CDenseVectorConstView  x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    CDenseVectorView       y(CArrayView(*N, Y, abs(*INCY)), *INCY<0);

    y = x;
}


void
BLAS(zcopy)(const INTEGER   *N,
            const cdouble   *X,
            const INTEGER   *INCX,
            cdouble         *Y,
            const INTEGER   *INCY)
{
    using std::abs;

    ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    ZDenseVectorView       y(ZArrayView(*N, Y, abs(*INCY)), *INCY<0);

    y = x;
}

} // extern "C"
