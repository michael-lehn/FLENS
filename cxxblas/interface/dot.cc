#include <cxxblas/cxxblas.cxx>

extern "C" {

#ifndef COMPLEX_FLOAT1
    typedef CBLAS_FLOAT                 CXXBLAS_FLOAT;
#else
    typedef std::complex<CBLAS_FLOAT>   CXXBLAS_FLOAT;
#endif

    CBLAS_FLOAT
    CBLAS_NAME(CBLAS_INT n,
               const CBLAS_FLOAT *_x, CBLAS_INT incX,
               const CBLAS_FLOAT *_y, CBLAS_INT incY)
#ifdef CREATE_CBLAS
    {
        const CXXBLAS_FLOAT *x = reinterpret_cast<const CXXBLAS_FLOAT *>(_x);
        const CXXBLAS_FLOAT *y = reinterpret_cast<const CXXBLAS_FLOAT *>(_y);

        CBLAS_FLOAT result;
        cxxblas::dot<CBLAS_INT,
                     CXXBLAS_FLOAT,
                     CBLAS_FLOAT>(n, x, incX, y, incY, result);
        return result;
    }
#else
    ;
#endif // CREATE_CBLAS

#ifdef CREATE_BLAS
    CBLAS_FLOAT
    BLAS_NAME(CBLAS_INT *_n,
              const CBLAS_FLOAT *x, CBLAS_INT *_incX,
              const CBLAS_FLOAT *y, CBLAS_INT *_incY)
    {

        CBLAS_INT n    = *_n;
        CBLAS_INT incX = *_incX;
        CBLAS_INT incY = *_incY;

        // the blas interface calls the cblas interface
        // so any blas-test will also test the cblas-interface
        return CBLAS_NAME(n, x, incX, y, incY);
    }
#endif // CREATE_BLAS

} // extern "C"