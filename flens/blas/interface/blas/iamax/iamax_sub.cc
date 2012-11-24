#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(isamax_sub)(const INTEGER   *N,
                 const float     *X,
                 const INTEGER   *INCX,
                 INTEGER         *RES)
{
    using std::abs;

    SDenseVectorConstView  x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    *RES = blas::iamax(x);
}

void
BLAS(idamax_sub)(const INTEGER   *N,
                 const double    *X,
                 const INTEGER   *INCX,
                 INTEGER         *RES)
{
    using std::abs;

    DDenseVectorConstView  x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    *RES = blas::iamax(x);
}

void
BLAS(icamax_sub)(const INTEGER   *N,
                 const cfloat    *X,
                 const INTEGER   *INCX,
                 INTEGER         *RES)
{
    using std::abs;

    CDenseVectorConstView  x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    *RES = blas::iamax(x);
}

void
BLAS(izamax_sub)(const INTEGER   *N,
                 const cdouble   *X,
                 const INTEGER   *INCX,
                 INTEGER         *RES)
{
    using std::abs;

    ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    *RES = blas::iamax(x);
}

} // extern "C"
