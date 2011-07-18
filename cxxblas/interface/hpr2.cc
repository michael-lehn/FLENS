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

#ifndef COMPLEX_FLOAT1
    typedef CBLAS_FLOAT                 CBLAS_ALPHA;
    typedef CBLAS_FLOAT                 CXXBLAS_FLOAT;
#else
    typedef const CBLAS_FLOAT *         CBLAS_ALPHA;
    typedef std::complex<CBLAS_FLOAT>   CXXBLAS_FLOAT;
#endif

    enum CBLAS_ORDER        {CblasRowMajor=101, CblasColMajor=102};
    enum CBLAS_UPLO         {CblasUpper=121, CblasLower=122};

    void
    CBLAS_NAME(enum CBLAS_ORDER _order,
               enum CBLAS_UPLO _upLo,
               CBLAS_INT n,
               CBLAS_ALPHA _alpha,
               const CBLAS_FLOAT *_x, CBLAS_INT incX,
               const CBLAS_FLOAT *_y, CBLAS_INT incY,
               CBLAS_FLOAT  *_A)
#ifdef CREATE_CBLAS
    {
        StorageOrder order = (_order==CblasColMajor) ? ColMajor
                                                     : RowMajor;
        StorageUpLo upLo = Lower;
        if (_upLo==CblasUpper) {
            upLo = Upper;
        }

        const CXXBLAS_FLOAT *x = reinterpret_cast<const CXXBLAS_FLOAT *>(_x);
        const CXXBLAS_FLOAT *y = reinterpret_cast<const CXXBLAS_FLOAT *>(_y);
        CXXBLAS_FLOAT *A = reinterpret_cast<CXXBLAS_FLOAT *>(_A);

    #ifndef COMPLEX_FLOAT1
        CXXBLAS_FLOAT alpha = _alpha;
    #else
        CXXBLAS_FLOAT alpha(_alpha[0], _alpha[1]);
    #endif

        CXXBLAS_FLOAT *__A = A;
        StorageOrder __order = order;

#   ifdef TEST_ROW_MAJOR
        __order = (order==ColMajor) ? RowMajor : ColMajor;
        allocatePackedStorage(n, __A);
        switchPackedStorageOrder(order, upLo, n, A, __A);
#   endif

        cxxblas::hpr2(__order, upLo,
                      n,
                      alpha,
                      x, incX,
                      y, incY,
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
              const CBLAS_FLOAT *y, const CBLAS_INT *_incY,
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
        CBLAS_INT incX    = *_incX;
        CBLAS_INT incY    = *_incY;

#   ifndef COMPLEX_FLOAT1
        CBLAS_ALPHA alpha  = *_alpha;
#   else
        CBLAS_ALPHA alpha = _alpha;
#   endif

        CBLAS_INT info = 0;
        if (incX==0) {
            info = 5;
        }
        if (incY==0) {
            info = 7;
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
                   y, incY,
                   A);
    }
#endif // CREATE_BLAS

} // extern "C"