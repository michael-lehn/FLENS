#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(srotmg)(float   *A,
             float   *B,
             float   *C,
             float   *D,
             float   *S)
{
#   ifdef TEST_DIRECT_CBLAS
        
        cblas_srotmg(A, B, C, D, S);
        
#   else
        SDenseVectorView s(SArrayView(5, S, 1), false);
        blas::rotmg(*A, *B, *C, *D, s);
        
#   endif
}

void
BLAS(drotmg)(double  *A,
             double  *B,
             double  *C,
             double  *D,
             double  *S)
{ 
#   ifdef TEST_DIRECT_CBLAS
        
        cblas_drotmg(A, B, C, D, S);
        
#   else
        DDenseVectorView s(DArrayView(5, S, 1), false);
        blas::rotmg(*A, *B, *C, *D, s);
        
#   endif
}

} // extern "C"
