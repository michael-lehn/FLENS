#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(saxpy)(const INTEGER   *N,
            const float     *ALPHA,
            const float     *X,
            const INTEGER   *INCX,
            float           *Y,
            const INTEGER   *INCY)
{
    using std::abs;

    SDenseVectorConstView  x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    SDenseVectorView       y(SArrayView(*N, Y, abs(*INCY)), *INCY<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    y += (*ALPHA) * x;
#   else
    blas::axpy(*ALPHA, x, y);
#   endif
}

void
BLAS(daxpy)(const INTEGER   *N,
            const double    *ALPHA,
            const double    *X,
            const INTEGER   *INCX,
            double          *Y,
            const INTEGER   *INCY)
{
    using std::abs;

    DDenseVectorConstView  x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    DDenseVectorView       y(DArrayView(*N, Y, abs(*INCY)), *INCY<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    y += (*ALPHA) * x;
#   else
    blas::axpy(*ALPHA, x, y);
#   endif
}

void
BLAS(caxpy)(const INTEGER   *N,
            const cfloat    *ALPHA,
            const cfloat    *X,
            const INTEGER   *INCX,
            cfloat          *Y,
            const INTEGER   *INCY)
{
    using std::abs;

    CDenseVectorConstView  x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    CDenseVectorView       y(CArrayView(*N, Y, abs(*INCY)), *INCY<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    y += (*ALPHA) * x;
#   else
    blas::axpy(*ALPHA, x, y);
#   endif
}

void
BLAS(zaxpy)(const INTEGER   *N,
            const cdouble   *ALPHA,
            const cdouble   *X,
            const INTEGER   *INCX,
            cdouble         *Y,
            const INTEGER   *INCY)
{
    using std::abs;

    ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    ZDenseVectorView       y(ZArrayView(*N, Y, abs(*INCY)), *INCY<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    y += (*ALPHA) * x;
#   else
    blas::axpy(*ALPHA, x, y);
#   endif
}

} // extern "C"
