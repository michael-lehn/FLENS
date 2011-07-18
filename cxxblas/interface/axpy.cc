#include <cxxblas/cxxblas.cxx>
#include <cstdlib>

extern "C" {

#ifndef COMPLEX_FLOAT1
    typedef CBLAS_FLOAT                 CBLAS_ALPHA;
    typedef CBLAS_FLOAT                 CXXBLAS_FLOAT;
#else
    typedef const CBLAS_FLOAT *         CBLAS_ALPHA;
    typedef std::complex<CBLAS_FLOAT>   CXXBLAS_FLOAT;
#endif

    void
    CBLAS_NAME(CBLAS_INT n,
               CBLAS_ALPHA _alpha, const CBLAS_FLOAT *_x, CBLAS_INT incX,
               CBLAS_FLOAT *_y, CBLAS_INT incY)
#ifdef CREATE_CBLAS
    {
        const CXXBLAS_FLOAT *x = reinterpret_cast<const CXXBLAS_FLOAT *>(_x);
        CXXBLAS_FLOAT *y = reinterpret_cast<CXXBLAS_FLOAT *>(_y);

#   ifndef COMPLEX_FLOAT1
        CXXBLAS_FLOAT alpha = _alpha;
#   else
        CXXBLAS_FLOAT alpha(_alpha[0], _alpha[1]);
#   endif

        cxxblas::axpy<CBLAS_INT,
                      CXXBLAS_FLOAT,
                      CXXBLAS_FLOAT>(n, alpha, x, incX, y, incY);
    }
#else
    ;
#endif // CREATE_CBLAS

#ifdef CREATE_BLAS
    void
    BLAS_NAME(CBLAS_INT *_n,
              CBLAS_FLOAT *_alpha, const CBLAS_FLOAT *x, CBLAS_INT *_incX,
              CBLAS_FLOAT *y, CBLAS_INT *_incY)
    {

        CBLAS_INT n       = *_n;
        CBLAS_INT incX    = *_incX;
        CBLAS_INT incY    = *_incY;
        
#   ifndef COMPLEX_FLOAT2
        CBLAS_ALPHA alpha  = *_alpha;
#   else
        CBLAS_ALPHA alpha = _alpha;
#   endif

        // the blas interface calls the cblas interface
        // so any blas-test will also test the cblas-interface
        CBLAS_NAME(n, alpha, x, incX, y, incY);
    }
#endif // CREATE_BLAS

} // extern "C"