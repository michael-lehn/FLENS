#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(sasum_sub)(const INTEGER   *N,
                const float     *X,
                const INTEGER   *INCX,
                float           *RES)
{
    using std::abs;

    SDenseVectorConstView  x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    blas::asum(x, *RES);
}


void
BLAS(dasum_sub)(const INTEGER   *N,
                const double    *X,
                const INTEGER   *INCX,
                double          *RES)
{
    using std::abs;

    DDenseVectorConstView  x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    blas::asum(x, *RES);
}

void
BLAS(scasum_sub)(const INTEGER   *N,
                 const cfloat    *X,
                 const INTEGER   *INCX,
                 float           *RES)
{
    using std::abs;

    CDenseVectorConstView  x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    blas::asum(x, *RES);
}


void
BLAS(dzasum_sub)(const INTEGER   *N,
                 const cdouble   *X,
                 const INTEGER   *INCX,
                 double          *RES)
{
    using std::abs;

    ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    blas::asum(x, *RES);
}


} // extern "C"
