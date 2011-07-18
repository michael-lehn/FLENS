#ifndef CXXLAPACK_INTERFACE_AUX_H
#define CXXLAPACK_INTERFACE_AUX_H 1

#include <complex>
#include <cxxlapack/cxxlapack.cxx>


#ifndef INT
#   define INT int
#endif


#if defined(SINGLE)
#   define FLOAT        float
#   define CXXFLOAT     float
#elif defined(DOUBLE)
#   define FLOAT        double
#   define CXXFLOAT     double
#elif defined(COMPLEX_SINGLE)
#   define FLOAT        float
#   define CXXFLOAT     std::complex<float>
#elif defined(COMPLEX_DOUBLE)
#   define FLOAT        double
#   define CXXFLOAT     std::complex<double>
#endif


using cxxblas::ColMajor;

using cxxblas::Transpose;
using cxxblas::NoTrans;
using cxxblas::Trans;
using cxxblas::ConjTrans;
using cxxblas::Conj;


extern "C" {

    #define XERBLA          xerbla_

    void
    XERBLA(const char *SRNAME, int *INFO, int len);

} // extern "C"

#endif // CXXLAPACK_INTERFACE_AUX_H