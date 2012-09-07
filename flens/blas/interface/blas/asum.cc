#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

float
BLAS(sasum)(const INTEGER   *N,
            const float     *X,
            const INTEGER   *INCX)
{
    using std::abs;

    SDenseVectorConstView  x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    double result;
    blas::asum(x, result);

    return result;
}


double
BLAS(dasum)(const INTEGER   *N,
            const double    *X,
            const INTEGER   *INCX)
{
    using std::abs;

    DDenseVectorConstView  x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    double result;
    blas::asum(x, result);

    return result;
}

float
BLAS(scasum)(const INTEGER   *N,
             const cfloat    *X,
             const INTEGER   *INCX)
{
    using std::abs;

    CDenseVectorConstView  x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    double result;
    blas::asum(x, result);

    return result;
}


double
BLAS(dzasum)(const INTEGER   *N,
             const cdouble   *X,
             const INTEGER   *INCX)
{
    using std::abs;

    ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    double result;
    blas::asum(x, result);

    return result;
}


} // extern "C"
