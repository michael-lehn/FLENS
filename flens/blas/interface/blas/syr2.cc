#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(ssyr2)(const char      *UPLO,
            const INTEGER   *N,
            const float     *ALPHA,
            const float     *X,
            const INTEGER   *INCX,
            const float     *Y,
            const INTEGER   *INCY,
            float           *_A,
            const INTEGER   *LDA)
{
    
#   ifdef TEST_DIRECT_CBLAS
        
        char    _UPLO   = toupper(*UPLO);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        
        cblas_ssyr2(CBLAS_ORDER::CblasColMajor,
                   cxxblas::CBLAS::getCblasType(upLo),
                   *N, *ALPHA,
                   X, *INCX, Y, *INCY, _A, *LDA);
        
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
            } else if (*INCY==0) {
                info = 7;
            } else if (*LDA<max(INTEGER(1),*N)) {
                info = 9;
            }
            if (info!=0) {
                BLAS(xerbla)("SSYR2 ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(_UPLO);

        SDenseVectorConstView x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        SDenseVectorConstView y(SConstArrayView(*N, Y, abs(*INCY)), *INCY<0);
        SSyMatrixView         A(SFullView(*N, *N, _A, *LDA), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;

            A += alpha*x*transpose(y) + alpha*y*transpose(x);
#       else
            blas::r2(*ALPHA, x, y, A);
#       endif
#   endif
}

void
BLAS(dsyr2)(const char      *UPLO,
            const INTEGER   *N,
            const double    *ALPHA,
            const double    *X,
            const INTEGER   *INCX,
            const double    *Y,
            const INTEGER   *INCY,
            double          *_A,
            const INTEGER   *LDA)
{
    
#   ifdef TEST_DIRECT_CBLAS
        
        char    _UPLO   = toupper(*UPLO);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        
        cblas_dsyr2(CBLAS_ORDER::CblasColMajor,
                   cxxblas::CBLAS::getCblasType(upLo),
                   *N, *ALPHA,
                   X, *INCX, Y, *INCY, _A, *LDA);
        
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
            } else if (*INCY==0) {
                info = 7;
            } else if (*LDA<max(INTEGER(1),*N)) {
                info = 9;
            }
            if (info!=0) {
                BLAS(xerbla)("DSYR2 ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(_UPLO);

        DDenseVectorConstView x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        DDenseVectorConstView y(DConstArrayView(*N, Y, abs(*INCY)), *INCY<0);
        DSyMatrixView         A(DFullView(*N, *N, _A, *LDA), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;

            A += alpha*x*transpose(y) + alpha*y*transpose(x);
#       else
            blas::r2(*ALPHA, x, y, A);
#       endif
#   endif
}

} // extern "C"
