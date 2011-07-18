#include <cxxblas/cxxblas.cxx>

extern "C" {

#ifndef COMPLEX_FLOAT1
    typedef CBLAS_FLOAT                 CXXBLAS_FLOAT;
#else
    typedef std::complex<CBLAS_FLOAT>   CXXBLAS_FLOAT;
#endif

#ifndef COMPLEX_FLOAT2
    typedef CBLAS_FLOAT                 CBLAS_ALPHA;
#else
    typedef const CBLAS_FLOAT *         CBLAS_ALPHA;
#endif

    void
    CBLAS_NAME(CBLAS_INT n, CBLAS_ALPHA _alpha,
               CBLAS_FLOAT *_x, CBLAS_INT incX)
#ifdef CREATE_CBLAS
    {
        CXXBLAS_FLOAT *x = reinterpret_cast<CXXBLAS_FLOAT *>(_x);


#   ifndef COMPLEX_FLOAT2
        CXXBLAS_FLOAT alpha = _alpha;
#   else
        CXXBLAS_FLOAT alpha(_alpha[0], _alpha[1]);
#   endif
        cxxblas::scal(n, alpha, x, incX);
    }
#else
    ;
#endif // CREATE_CBLAS

#ifdef CREATE_BLAS
    void
    BLAS_NAME(const CBLAS_INT *_n, const CBLAS_FLOAT *_alpha,
              CBLAS_FLOAT *x, const CBLAS_INT *_incX)
    {

        CBLAS_INT n         = *_n;
        CBLAS_INT incX      = *_incX;

#   ifndef COMPLEX_FLOAT2
        CBLAS_ALPHA alpha  = *_alpha;
#   else
        CBLAS_ALPHA alpha = _alpha;
#   endif

        // the blas interface calls the cblas interface
        // so any blas-test will also test the cblas-interface
        CBLAS_NAME(n, alpha, x, incX);
    }
#endif // CREATE_BLAS

} // extern "C"