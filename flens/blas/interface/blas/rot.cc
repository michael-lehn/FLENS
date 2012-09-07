#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(srot)(const INTEGER   *N,
           float           *X,
           const INTEGER   *INCX,
           float           *Y,
           const INTEGER   *INCY,
           const float     *C,
           const float     *S)
{
    using std::abs;

    SDenseVectorView  x(SArrayView(*N, X, abs(*INCX)), *INCX<0);
    SDenseVectorView  y(SArrayView(*N, Y, abs(*INCY)), *INCY<0);

    blas::rot(x, y, *C, *S);
}

void
BLAS(drot)(const INTEGER   *N,
           double          *X,
           const INTEGER   *INCX,
           double          *Y,
           const INTEGER   *INCY,
           const double    *C,
           const double    *S)
{
    using std::abs;

    DDenseVectorView  x(DArrayView(*N, X, abs(*INCX)), *INCX<0);
    DDenseVectorView  y(DArrayView(*N, Y, abs(*INCY)), *INCY<0);

    blas::rot(x, y, *C, *S);
}

void
BLAS(crot)(const INTEGER   *N,
           cfloat          *X,
           const INTEGER   *INCX,
           cfloat          *Y,
           const INTEGER   *INCY,
           const cfloat    *C,
           const cfloat    *S)
{
    using std::abs;

    CDenseVectorView  x(CArrayView(*N, X, abs(*INCX)), *INCX<0);
    CDenseVectorView  y(CArrayView(*N, Y, abs(*INCY)), *INCY<0);

    blas::rot(x, y, *C, *S);
}

void
BLAS(zrot)(const INTEGER   *N,
           cdouble         *X,
           const INTEGER   *INCX,
           cdouble         *Y,
           const INTEGER   *INCY,
           const cdouble   *C,
           const cdouble   *S)
{
    using std::abs;

    ZDenseVectorView  x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);
    ZDenseVectorView  y(ZArrayView(*N, Y, abs(*INCY)), *INCY<0);

    blas::rot(x, y, *C, *S);
}


} // extern "C"
