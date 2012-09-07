#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(sscal)(const INTEGER   *N,
            const float     *ALPHA,
            float           *X,
            const INTEGER   *INCX)
{
    using std::abs;

    SDenseVectorView  x(SArrayView(*N, X, abs(*INCX)), *INCX<0);

    x *= (*ALPHA);
}

void
BLAS(dscal)(const INTEGER   *N,
            const double    *ALPHA,
            double          *X,
            const INTEGER   *INCX)
{
    using std::abs;

    DDenseVectorView  x(DArrayView(*N, X, abs(*INCX)), *INCX<0);

    x *= (*ALPHA);
}

void
BLAS(cscal)(const INTEGER   *N,
            const cfloat    *ALPHA,
            cfloat          *X,
            const INTEGER   *INCX)
{
    using std::abs;

    CDenseVectorView  x(CArrayView(*N, X, abs(*INCX)), *INCX<0);

    x *= (*ALPHA);
}

void
BLAS(csscal)(const INTEGER   *N,
             const float     *ALPHA,
             cfloat          *X,
             const INTEGER   *INCX)
{
    using std::abs;

    CDenseVectorView  x(CArrayView(*N, X, abs(*INCX)), *INCX<0);

    x *= (*ALPHA);
}

void
BLAS(zscal)(const INTEGER   *N,
            const cdouble   *ALPHA,
            cdouble         *X,
            const INTEGER   *INCX)
{
    using std::abs;

    ZDenseVectorView  x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);

    x *= (*ALPHA);
}

void
BLAS(zdscal)(const INTEGER   *N,
             const double    *ALPHA,
             cdouble         *X,
             const INTEGER   *INCX)
{
    using std::abs;

    ZDenseVectorView  x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);

    x *= (*ALPHA);
}

} // extern "C"
