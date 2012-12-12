#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(srotm)(const INTEGER   *N,
            float           *X,
            const INTEGER   *INCX,
            float           *Y,
            const INTEGER   *INCY,
            const float     *P)
{
#   ifdef TEST_DIRECT_CBLAS
    
        cblas_srotm(*N, X, *INCX, Y, *INCY, *P);
    
#   else
        using std::abs;
    
        SDenseVectorView  x(SArrayView(*N, X, abs(*INCX)), *INCX<0);
        SDenseVectorView  y(SArrayView(*N, Y, abs(*INCY)), *INCY<0);
        SDenseVectorConstView p(SConstArrayView(5, P, 1), false);
    
        blas::rotm(x, y, p);
        
#   endif
    
}

void
BLAS(drotm)(const INTEGER   *N,
            double          *X,
            const INTEGER   *INCX,
            double          *Y,
            const INTEGER   *INCY,
            const double    *P)
{
#   ifdef TEST_DIRECT_CBLAS
        
        cblas_drotm(*N, X, *INCX, Y, *INCY, *P);
        
#   else
        
        DDenseVectorView  x(DArrayView(*N, X, abs(*INCX)), *INCX<0);
        DDenseVectorView  y(DArrayView(*N, Y, abs(*INCY)), *INCY<0);
        DDenseVectorConstView p(DConstArrayView(5, P, 1), false);
    
        blas::rotm(x, y, p);
    
#   endif
}

} // extern "C"
