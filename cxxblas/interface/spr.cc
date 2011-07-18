#include <cxxblas/cxxblas.cxx>
#include <cxxblas/interface/aux.h>

using cxxblas::StorageOrder;
using cxxblas::ColMajor;
using cxxblas::RowMajor;

using cxxblas::StorageUpLo;
using cxxblas::Upper;
using cxxblas::Lower;

extern "C" {

    void xerbla_(const char* srname, int* info);

    enum CBLAS_ORDER        {CblasRowMajor=101, CblasColMajor=102};
    enum CBLAS_UPLO         {CblasUpper=121, CblasLower=122};

    void
    CBLAS_NAME(enum CBLAS_ORDER _order,
               enum CBLAS_UPLO _upLo,
               CBLAS_INT n,
               CBLAS_FLOAT alpha,
               const CBLAS_FLOAT  *x, CBLAS_INT incX,
               CBLAS_FLOAT  *A)
#ifdef CREATE_CBLAS
    {
        StorageOrder order = (_order==CblasColMajor) ? ColMajor
                                                     : RowMajor;
        StorageUpLo upLo = Lower;
        if (_upLo==CblasUpper) {
            upLo = Upper;
        }

        CBLAS_FLOAT *__A = A;
        StorageOrder __order = order;

#   ifdef TEST_ROW_MAJOR
        __order = (order==ColMajor) ? RowMajor : ColMajor;
        allocatePackedStorage(n, __A);
        switchPackedStorageOrder(order, upLo, n, A, __A);
#   endif

        cxxblas::spr(__order, upLo,
                     n,
                     alpha,
                     x, incX,
                     __A);

#   ifdef TEST_ROW_MAJOR
        switchPackedStorageOrder(__order, upLo, n, __A, A);
        releaseStorage(__A);
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
              const CBLAS_FLOAT *x, const CBLAS_INT *_incX,
              CBLAS_FLOAT *A)
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
        CBLAS_INT incX    = *_incX;

        CBLAS_INT info = 0;
        if (incX==0) {
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
        CBLAS_NAME(CblasColMajor, upLo,
                   n,
                   alpha,
                   x, incX,
                   A);
    }
#endif // CREATE_BLAS

} // extern "C"