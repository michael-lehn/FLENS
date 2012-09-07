#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

INTEGER
BLAS(isamax)(const INTEGER   *N,
             const float     *X,
             const INTEGER   *INCX)
{
    using std::abs;

    SDenseVectorConstView  x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    return blas::iamax(x);
}

INTEGER
BLAS(idamax)(const INTEGER   *N,
             const double    *X,
             const INTEGER   *INCX)
{
    using std::abs;

    DDenseVectorConstView  x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    return blas::iamax(x);
}

INTEGER
BLAS(icamax)(const INTEGER   *N,
             const cfloat    *X,
             const INTEGER   *INCX)
{
    using std::abs;

    CDenseVectorConstView  x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    return blas::iamax(x);
}

INTEGER
BLAS(izamax)(const INTEGER   *N,
             const cdouble   *X,
             const INTEGER   *INCX)
{
    using std::abs;

    ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    return blas::iamax(x);
}

} // extern "C"
