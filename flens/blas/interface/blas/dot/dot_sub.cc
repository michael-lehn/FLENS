#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(sdot_sub)(const INTEGER   *N,
               const float     *X,
               const INTEGER   *INCX,
               const float     *Y,
               const INTEGER   *INCY,
               float           *RES)
{
#   ifdef TEST_DIRECT_CBLAS

        *RES = cblas_sdot(*N, X, *INCX, Y, *INCY);

#   else

        using std::abs;

        SDenseVectorConstView x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        SDenseVectorConstView y(SConstArrayView(*N, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            *RES = x*y;
#       else
            *RES = blas::dot(x, y);
#       endif
#   endif
}

void
BLAS(ddot_sub)(const INTEGER   *N,
               const double    *X,
               const INTEGER   *INCX,
               const double    *Y,
               const INTEGER   *INCY,
               double          *RES)
{
#   ifdef TEST_DIRECT_CBLAS

        *RES = cblas_ddot(*N, X, *INCX, Y, *INCY);

#   else
        using std::abs;

        DDenseVectorConstView  x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        DDenseVectorConstView  y(DConstArrayView(*N, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            *RES = x*y;
#       else
            *RES = blas::dot(x, y);
#       endif
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
#   ifdef TEST_DIRECT_CBLAS

    cblas_cdotc_sub(*N,
                    reinterpret_cast<const float *>(X), *INCX,
                    reinterpret_cast<const float *>(Y), *INCY,
                    reinterpret_cast<float *>(RES));

#   else
        using std::abs;

        CDenseVectorConstView x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        CDenseVectorConstView y(CConstArrayView(*N, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            *RES = conjugate(x)*y;
#       else
            *RES = blas::dot(x, y);
#       endif
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
#   ifdef TEST_DIRECT_CBLAS

    cblas_cdotu_sub(*N,
                    reinterpret_cast<const float *>(X), *INCX,
                    reinterpret_cast<const float *>(Y), *INCY,
                    reinterpret_cast<float *>(RES));

#   else
        using std::abs;

        CDenseVectorConstView x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        CDenseVectorConstView y(CConstArrayView(*N, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            *RES = x*y;
#       else
            *RES = blas::dotu(x, y);
#       endif
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
#   ifdef TEST_DIRECT_CBLAS

    cblas_zdotc_sub(*N,
                    reinterpret_cast<const double *>(X), *INCX,
                    reinterpret_cast<const double *>(Y), *INCY,
                    reinterpret_cast<double *>(RES));

#   else
        using std::abs;

        ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        ZDenseVectorConstView  y(ZConstArrayView(*N, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            *RES = conjugate(x)*y;
#       else
            *RES = blas::dot(x, y);
#       endif
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

#   ifdef TEST_DIRECT_CBLAS

    cblas_zdotu_sub(*N,
                    reinterpret_cast<const double *>(X), *INCX,
                    reinterpret_cast<const double *>(Y), *INCY,
                    reinterpret_cast<double *>(RES));

#   else
        using std::abs;

        ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        ZDenseVectorConstView  y(ZConstArrayView(*N, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            *RES = x*y;
#       else
            *RES = blas::dotu(x, y);
#       endif
#   endif

}

void
BLAS(sdsdot_sub)(const INTEGER   *N,
                 const float     *B,
                 const float     *X,
                 const INTEGER   *INCX,
                 const float     *Y,
                 const INTEGER   *INCY,
                 float           *RES)
{
#   ifdef TEST_DIRECT_CBLAS

        *RES = cblas_sdsdot(*N, *B, X, *INCX, Y, *INCY);

#   else
        using std::abs;

        SDenseVectorConstView x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        SDenseVectorConstView y(SConstArrayView(*N, Y, abs(*INCY)), *INCY<0);

        double tmp;
        blas::dot(x, y, tmp);

        tmp += double(*B);
        *RES = tmp;
#   endif
}

void
BLAS(dsdot_sub)(const INTEGER   *N,
                const float     *X,
                const INTEGER   *INCX,
                const float     *Y,
                const INTEGER   *INCY,
                double          *RES)
{
#   ifdef TEST_DIRECT_CBLAS

    *RES =  cblas_dsdot(*N, X, *INCX, Y, *INCY);

#   else

        using std::abs;

        SDenseVectorConstView x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        SDenseVectorConstView y(SConstArrayView(*N, Y, abs(*INCY)), *INCY<0);

        blas::dot(x, y, *RES);

#   endif
}


} // extern "C"
