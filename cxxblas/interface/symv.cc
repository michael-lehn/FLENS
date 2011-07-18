#include <cxxblas/cxxblas.cxx>
#include <cxxblas/interface/aux.h>

using cxxblas::StorageOrder;
using cxxblas::ColMajor;
using cxxblas::RowMajor;
using cxxblas::Transpose;
using cxxblas::NoTrans;
using cxxblas::Trans;
using cxxblas::Conj;
using cxxblas::ConjTrans;
using cxxblas::StorageUpLo;
using cxxblas::Upper;
using cxxblas::Lower;


extern "C" {

    void xerbla_(const char* srname, int* info);

    enum CBLAS_ORDER        {CblasRowMajor=101, CblasColMajor=102};
    enum CBLAS_TRANSPOSE    {CblasNoTrans=111, CblasTrans=112,
                             CblasConjTrans=113, CblasConjNoTrans=114};
    enum CBLAS_UPLO         {CblasUpper=121, CblasLower=122};

    void
    CBLAS_NAME(enum CBLAS_ORDER _order,
               enum CBLAS_UPLO _upLo,
               CBLAS_INT n,
               CBLAS_FLOAT alpha,
               const CBLAS_FLOAT  *A, CBLAS_INT ldA,
               const CBLAS_FLOAT  *x, CBLAS_INT incX,
               CBLAS_FLOAT beta,
               CBLAS_FLOAT  *y, CBLAS_INT incY)
#ifdef CREATE_CBLAS
    {
        StorageOrder order = (_order==CblasColMajor) ? ColMajor
                                                     : RowMajor;
        StorageUpLo upLo = Lower;
        if (_upLo==CblasUpper) {
            upLo = Upper;
        }

#   ifdef TEST_ROW_MAJOR
        switchFullStorageOrder(order, n, n, A, ldA);
        order = (order==ColMajor) ? RowMajor : ColMajor;
#   endif

        cxxblas::symv(order,
                      upLo,
                      n,
                      alpha,
                      A, ldA,
                      x, incX,
                      beta,
                      y, incY);

#   ifdef TEST_ROW_MAJOR
        releaseStorage(A);
#   endif
    }
#else
    ;
#endif // CREATE_CBLAS

#ifdef CREATE_BLAS
    void
    BLAS_NAME(const char *_upLo,
              const CBLAS_INT *_n,
              const CBLAS_FLOAT *_alpha,
              const CBLAS_FLOAT *A, const CBLAS_INT *_ldA,
              const CBLAS_FLOAT *x, const CBLAS_INT *_incX,
              const CBLAS_FLOAT *_beta,
              CBLAS_FLOAT *y, const CBLAS_INT *_incY)
    {
        bool checkUpLo = false;
        CBLAS_UPLO upLo = CblasUpper;
        if ((*_upLo=='L') || (*_upLo=='l')) {
            upLo = CblasLower;
            checkUpLo = true;
        }
        if ((*_upLo=='U') || (*_upLo=='u')) {
            checkUpLo = true;
        }

        CBLAS_INT n       = *_n;
        CBLAS_FLOAT alpha = *_alpha;
        CBLAS_INT ldA     = *_ldA;
        CBLAS_INT incX    = *_incX;
        CBLAS_FLOAT beta  = *_beta;
        CBLAS_INT incY    = *_incY;

        CBLAS_INT info = 0;
        if (incY==0) {
            info = 10;
        }
        if (incX==0) {
            info = 7; 
        }
        if (ldA<std::max(1, n)) {
            info = 5;
        }
        if (n<0) {
            info = 2;
        }
        if (!checkUpLo) {
            info = 1;
        }
        if (info!=0) {
            char blasName[6];
            strncpy(blasName, BLAS_NAME_STR, 6);
            for (int i=0; i<6; ++i) {
                blasName[i] = std::toupper(blasName[i]);
            }
            xerbla_(blasName, &info);
          return;
        }

        // the blas interface calls the cblas interface
        // so any blas-test will also test the cblas-interface
        CBLAS_NAME(CblasColMajor,
                   upLo,
                   n,
                   alpha,
                   A, ldA,
                   x, incX,
                   beta,
                   y, incY);
    }
#endif // CREATE_BLAS

} // extern "C"