#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(sdot_sub)(const INTEGER   *N,
               const float    *X,
               const INTEGER  *INCX,
               const float    *Y,
               const INTEGER  *INCY,
               float          *RES)
{
    using std::abs;

    SDenseVectorConstView x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    SDenseVectorConstView y(SConstArrayView(*N, Y, abs(*INCY)), *INCY<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    *RES = x*y;
#   else
    *RES = blas::dot(x, y);
#   endif
}

void
BLAS(ddot_sub)(const INTEGER   *N,
               const double   *X,
               const INTEGER  *INCX,
               const double   *Y,
               const INTEGER  *INCY,
               double         *RES)
{
    using std::abs;

    DDenseVectorConstView  x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    DDenseVectorConstView  y(DConstArrayView(*N, Y, abs(*INCY)), *INCY<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    *RES = x*y;
#   else
    *RES = blas::dot(x, y);
#   endif
}

void
BLAS(cdotc_sub)(const INTEGER  *N,
                const cfloat   *X,
                const INTEGER  *INCX,
                const cfloat   *Y,
                const INTEGER  *INCY,
                cfloat         *RES)
{
    using std::abs;

    CDenseVectorConstView x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    CDenseVectorConstView y(CConstArrayView(*N, Y, abs(*INCY)), *INCY<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    *RES = conjugate(x)*y;
#   else
    *RES = blas::dot(x, y);
#   endif
}

void
BLAS(cdotu_sub)(const INTEGER  *N,
                const cfloat   *X,
                const INTEGER  *INCX,
                const cfloat   *Y,
                const INTEGER  *INCY,
                cfloat         *RES)
{
    using std::abs;

    CDenseVectorConstView x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    CDenseVectorConstView y(CConstArrayView(*N, Y, abs(*INCY)), *INCY<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    *RES = x*y;
#   else
    *RES = blas::dotu(x, y);
#   endif
}

void
BLAS(zdotc_sub)(const INTEGER  *N,
                const cdouble  *X,
                const INTEGER  *INCX,
                const cdouble  *Y,
                const INTEGER  *INCY,
                cdouble        *RES)
{
    using std::abs;

    ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    ZDenseVectorConstView  y(ZConstArrayView(*N, Y, abs(*INCY)), *INCY<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    *RES = conjugate(x)*y;
#   else
    *RES = blas::dot(x, y);
#   endif
}

void
BLAS(zdotu_sub)(const INTEGER  *N,
                const cdouble  *X,
                const INTEGER  *INCX,
                const cdouble  *Y,
                const INTEGER  *INCY,
                cdouble        *RES)
{
    using std::abs;

    ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    ZDenseVectorConstView  y(ZConstArrayView(*N, Y, abs(*INCY)), *INCY<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    *RES = x*y;
#   else
    *RES = blas::dotu(x, y);
#   endif
}

} // extern "C"
