#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(cher)(const char      *UPLO,
           const INTEGER   *N,
           const float     *ALPHA,
           const cfloat    *X,
           const INTEGER   *INCX,
           cfloat          *_A,
           const INTEGER   *LDA)
{
    
#   ifdef TEST_DIRECT_CBLAS
    
        char    _UPLO   = toupper(*UPLO);
    
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
    
        cblas_cher(CBLAS_ORDER::CblasColMajor,
                   cxxblas::CBLAS::getCblasType(upLo),
                   *N,
                   *ALPHA,
                   reinterpret_cast<const float *>(X), *INCX,
                   reinterpret_cast<float *>(_A), *LDA);
    
#   else
    
        using std::abs;
        using std::max;
    
        char    _UPLO = toupper(*UPLO);

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;

            if (_UPLO!='U' && _UPLO!='L') {
                info = 1;
            } else if (*N<0) {
                info = 2;
            } else if (*INCX==0) {
                info = 5;
            } else if (*LDA<max(INTEGER(1),*N)) {
                info = 7;
            }
            if (info!=0) {
                BLAS(xerbla)("CHER  ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(_UPLO);

        CDenseVectorConstView x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        CHeMatrixView         A(CFullView(*N, *N, _A, *LDA), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            A += alpha*x*conjTrans(x);
#       else
            blas::r(*ALPHA, x, A);
#       endif
#   endif
}

void
BLAS(zher)(const char      *UPLO,
           const INTEGER   *N,
           const double    *ALPHA,
           const cdouble   *X,
           const INTEGER   *INCX,
           cdouble         *_A,
           const INTEGER   *LDA)
{
#   ifdef TEST_DIRECT_CBLAS
    
        char    _UPLO   = toupper(*UPLO);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        
        cblas_zher(CBLAS_ORDER::CblasColMajor,
                   cxxblas::CBLAS::getCblasType(upLo),
                   *N,
                   *ALPHA,
                   reinterpret_cast<const double *>(X), *INCX,
                   reinterpret_cast<double *>(_A), *LDA);
        
#   else
    
        using std::abs;
        using std::max;

        char    _UPLO = toupper(*UPLO);
    
#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;

            if (_UPLO!='U' && _UPLO!='L') {
                info = 1;
            } else if (*N<0) {
                info = 2;
            } else if (*INCX==0) {
                info = 5;
            } else if (*LDA<max(INTEGER(1),*N)) {
                info = 7;
            }
            if (info!=0) {
                BLAS(xerbla)("ZHER  ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(_UPLO);

        ZDenseVectorConstView x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        ZHeMatrixView         A(ZFullView(*N, *N, _A, *LDA), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            A += alpha*x*conjTrans(x);
#       else
            blas::r(*ALPHA, x, A);
#       endif
#   endif
}

} // extern "C"
