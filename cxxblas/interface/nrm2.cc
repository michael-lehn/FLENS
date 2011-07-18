#include <cxxblas/cxxblas.cxx>

extern "C" {

#ifndef COMPLEX_FLOAT1
    typedef CBLAS_FLOAT                 CXXBLAS_FLOAT;
#else
    typedef std::complex<CBLAS_FLOAT>   CXXBLAS_FLOAT;
#endif

    CBLAS_FLOAT
    CBLAS_NAME(CBLAS_INT n, const CBLAS_FLOAT *_x, CBLAS_INT incX)
#ifdef CREATE_CBLAS
    {
        const CXXBLAS_FLOAT *x = reinterpret_cast<const CXXBLAS_FLOAT *>(_x);

        CBLAS_FLOAT result;
        cxxblas::nrm2<CBLAS_INT,
                      CXXBLAS_FLOAT,
                      CBLAS_FLOAT>(n, x, incX, result);
        return result;
    }
#else
    ;
#endif // CREATE_CBLAS

#ifdef CREATE_BLAS
    CBLAS_FLOAT
    BLAS_NAME(CBLAS_INT *_n, const CBLAS_FLOAT *x, CBLAS_INT *_incX)
    {

        CBLAS_INT n    = *_n;
        CBLAS_INT incX = *_incX;

        // the blas interface calls the cblas interface
        // so any blas-test will also test the cblas-interface
        return CBLAS_NAME(n, x, incX);
    }
#endif // CREATE_BLAS

} // extern "C"