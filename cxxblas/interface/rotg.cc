#include <cxxblas/cxxblas.cxx>

extern "C" {

#ifndef COMPLEX_FLOAT1
    typedef CBLAS_FLOAT                 CXXBLAS_FLOAT;
#else
    typedef std::complex<CBLAS_FLOAT>   CXXBLAS_FLOAT;
#endif

    void
    CBLAS_NAME(CBLAS_FLOAT *a, CBLAS_FLOAT *b,
               CBLAS_FLOAT *c, CBLAS_FLOAT *s)
#ifdef CREATE_CBLAS
    {
        cxxblas::rotg<CXXBLAS_FLOAT,
                      CXXBLAS_FLOAT,
                      CXXBLAS_FLOAT>(*a, *b, *c, *s);
    }
#else
    ;
#endif // CREATE_CBLAS

#ifdef CREATE_BLAS
    void
    BLAS_NAME(CBLAS_FLOAT *a, CBLAS_FLOAT *b,
              CBLAS_FLOAT *c, CBLAS_FLOAT *s)
    {
        // the blas interface calls the cblas interface
        // so any blas-test will also test the cblas-interface
        CBLAS_NAME(a, b, c, s);
    }
#endif // CREATE_BLAS

} // extern "C"