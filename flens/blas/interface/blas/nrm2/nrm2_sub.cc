#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(snrm2_sub)(const INTEGER   *N,
                const float     *X,
                const INTEGER   *INCX,
                float           *RES)
{
    using std::abs;

    SDenseVectorConstView  x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    *RES = blas::nrm2(x);
}

void
BLAS(dnrm2_sub)(const INTEGER   *N,
                const double    *X,
                const INTEGER   *INCX,
                double          *RES)
{
    using std::abs;

    DDenseVectorConstView  x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    *RES = blas::nrm2(x);
}

void
BLAS(scnrm2_sub)(const INTEGER   *N,
                 const cfloat    *X,
                 const INTEGER   *INCX,
                 float           *RES)
{
    using std::abs;

    CDenseVectorConstView  x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    *RES = blas::nrm2(x);
}

void
BLAS(dznrm2_sub)(const INTEGER   *N,
                 const cdouble   *X,
                 const INTEGER   *INCX,
                 double          *RES)
{
    using std::abs;

    ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    *RES = blas::nrm2(x);
}

} // extern "C"
