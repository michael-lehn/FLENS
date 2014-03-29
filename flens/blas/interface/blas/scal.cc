#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(sscal)(const INTEGER   *N,
            const float     *ALPHA,
            float           *X,
            const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS

        cblas_sscal(*N, *ALPHA, X, *INCX);

#   else

        using std::abs;

        SDenseVectorView  x(SArrayView(*N, X, abs(*INCX)), *INCX<0);

#       ifdef TEST_OVERLOADED_OPERATORS

            x *= (*ALPHA);

#       else

            blas::scal(*ALPHA, x);

#       endif

#   endif
}

void
BLAS(dscal)(const INTEGER   *N,
            const double    *ALPHA,
            double          *X,
            const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS

        cblas_dscal(*N, *ALPHA, X, *INCX);

#   else

        using std::abs;

        DDenseVectorView  x(DArrayView(*N, X, abs(*INCX)), *INCX<0);

#       ifdef TEST_OVERLOADED_OPERATORS

            x *= (*ALPHA);

#       else

            blas::scal(*ALPHA, x);

#       endif
#   endif
}

void
BLAS(cscal)(const INTEGER   *N,
            const cfloat    *ALPHA,
            cfloat          *X,
            const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS

        cblas_cscal(*N,
                    reinterpret_cast<const float *>(ALPHA),
                    reinterpret_cast<float *>(X), *INCX);

#   else
        using std::abs;

        CDenseVectorView  x(CArrayView(*N, X, abs(*INCX)), *INCX<0);

#       ifdef TEST_OVERLOADED_OPERATORS

            x *= (*ALPHA);

#       else

            blas::scal(*ALPHA, x);

#       endif

#   endif
}

void
BLAS(csscal)(const INTEGER   *N,
             const float     *ALPHA,
             cfloat          *X,
             const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS

    cblas_csscal(*N,
                *ALPHA,
                reinterpret_cast<float *>(X), *INCX);

#   else

        using std::abs;

        CDenseVectorView  x(CArrayView(*N, X, abs(*INCX)), *INCX<0);

#       ifdef TEST_OVERLOADED_OPERATORS

            x *= (*ALPHA);

#       else

            blas::scal(*ALPHA, x);

#       endif
#   endif
}

void
BLAS(zscal)(const INTEGER   *N,
            const cdouble   *ALPHA,
            cdouble         *X,
            const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS

    cblas_zscal(*N,
                reinterpret_cast<const double *>(ALPHA),
                reinterpret_cast<double *>(X), *INCX);

#   else

        using std::abs;

        ZDenseVectorView  x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);

#       ifdef TEST_OVERLOADED_OPERATORS

            x *= (*ALPHA);

#       else

            blas::scal(*ALPHA, x);

#       endif

#   endif
}

void
BLAS(zdscal)(const INTEGER   *N,
             const double    *ALPHA,
             cdouble         *X,
             const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS

    cblas_cscal(*N,
                reinterpret_cast<const float *>(ALPHA),
                reinterpret_cast<float *>(X), *INCX);

#   else

        using std::abs;

        ZDenseVectorView  x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);

#       ifdef TEST_OVERLOADED_OPERATORS

            x *= (*ALPHA);

#       else

            blas::scal(*ALPHA, x);

#       endif
#   endif
}

} // extern "C"
