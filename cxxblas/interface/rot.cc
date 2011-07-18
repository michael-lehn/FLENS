#include <cxxblas/cxxblas.cxx>

extern "C" {

#ifndef COMPLEX_FLOAT1
    typedef CBLAS_FLOAT                 CXXBLAS_FLOAT;
#else
    typedef std::complex<CBLAS_FLOAT>   CXXBLAS_FLOAT;
#endif

    void
    CBLAS_NAME(CBLAS_INT n,
               CBLAS_FLOAT *_x, CBLAS_INT incX,
               CBLAS_FLOAT *_y, CBLAS_INT incY,
               CBLAS_FLOAT c, CBLAS_FLOAT s)
#ifdef CREATE_CBLAS
    {
        CXXBLAS_FLOAT *x = reinterpret_cast<CXXBLAS_FLOAT *>(_x);
        CXXBLAS_FLOAT *y = reinterpret_cast<CXXBLAS_FLOAT *>(_y);

        cxxblas::rot<CBLAS_INT,
                     CXXBLAS_FLOAT,
                     CXXBLAS_FLOAT,
                     CXXBLAS_FLOAT>(n, x, incX, y, incY, c, s);
    }
#else
    ;
#endif // CREATE_CBLAS

#ifdef CREATE_BLAS
    void
    BLAS_NAME(CBLAS_INT *_n,
              CBLAS_FLOAT *x, CBLAS_INT *_incX,
              CBLAS_FLOAT *y, CBLAS_INT *_incY,
              CBLAS_FLOAT *_c, CBLAS_FLOAT *_s)
    {

        CBLAS_INT n    = *_n;
        CBLAS_INT incX = *_incX;
        CBLAS_INT incY = *_incY;
        CBLAS_FLOAT c  = *_c; 
        CBLAS_FLOAT s  = *_s; 

        // the blas interface calls the cblas interface
        // so any blas-test will also test the cblas-interface
        CBLAS_NAME(n, x, incX, y, incY, c, s);
    }
#endif // CREATE_BLAS

} // extern "C"