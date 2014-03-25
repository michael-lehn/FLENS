#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(sasum_sub)(const INTEGER   *N,
                const float     *X,
                const INTEGER   *INCX,
                float           *RES)
{
#   ifdef TEST_DIRECT_CBLAS
    
        (*RES) = cblas_sasum(*N, X, *INCX);
    
#   else
        using std::abs;

        SDenseVectorConstView  x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);

        blas::asum(x, *RES);
#   endif
}


void
BLAS(dasum_sub)(const INTEGER   *N,
                const double    *X,
                const INTEGER   *INCX,
                double          *RES)
{
#   ifdef TEST_DIRECT_CBLAS
        
        (*RES) = cblas_dasum(*N, X, *INCX);
        
#   else
        using std::abs;

        DDenseVectorConstView  x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);

        blas::asum(x, *RES);
#   endif
}

void
BLAS(scasum_sub)(const INTEGER   *N,
                 const cfloat    *X,
                 const INTEGER   *INCX,
                 float           *RES)
{
#   ifdef TEST_DIRECT_CBLAS
    
        (*RES) = cblas_scasum(*N, reinterpret_cast<const float *>(X), *INCX);
    
#   else
    
        using std::abs;

        CDenseVectorConstView  x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);

        blas::asum(x, *RES);
    
#   endif
}


void
BLAS(dzasum_sub)(const INTEGER   *N,
                 const cdouble   *X,
                 const INTEGER   *INCX,
                 double          *RES)
{
#   ifdef TEST_DIRECT_CBLAS
        
        (*RES) = cblas_dzasum(*N, reinterpret_cast<const double *>(X), *INCX);
        
#   else
        using std::abs;

        ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);

        blas::asum(x, *RES);
    
#   endif
    
}


} // extern "C"
