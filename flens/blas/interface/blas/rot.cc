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
#   ifdef TEST_DIRECT_CBLAS

        cblas_srot(*N, X, *INCX, Y, *INCY, *C, *S);

#   else

        using std::abs;

        SDenseVectorView  x(SArrayView(*N, X, abs(*INCX)), *INCX<0);
        SDenseVectorView  y(SArrayView(*N, Y, abs(*INCY)), *INCY<0);

        blas::rot(x, y, *C, *S);

#   endif
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
#   ifdef TEST_DIRECT_CBLAS

        cblas_drot(*N, X, *INCX, Y, *INCY, *C, *S);

#   else

        using std::abs;

        DDenseVectorView  x(DArrayView(*N, X, abs(*INCX)), *INCX<0);
        DDenseVectorView  y(DArrayView(*N, Y, abs(*INCY)), *INCY<0);

        blas::rot(x, y, *C, *S);

#   endif

}

//
//  crot and zrot are the LAPACK auxiliary routines CROT and ZROT: a plane
//  rotation with real cosine C and complex sine S,
//
//      x <- C*x + S*y,    y <- C*y - conj(S)*x.
//
void
BLAS(crot)(const INTEGER   *N,
           cfloat          *X,
           const INTEGER   *INCX,
           cfloat          *Y,
           const INTEGER   *INCY,
           const float     *C,
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
           const double    *C,
           const cdouble   *S)
{
    using std::abs;

    ZDenseVectorView  x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);
    ZDenseVectorView  y(ZArrayView(*N, Y, abs(*INCY)), *INCY<0);

    blas::rot(x, y, *C, *S);
}

} // extern "C"
