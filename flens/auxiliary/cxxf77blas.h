#ifndef FLENS_AUXILIARY_CXXF77BLAS_H
#define FLENS_AUXILIARY_CXXF77BLAS_H 1

#include <flens/typedefs.h>
#include <type_traits>

namespace flens { namespace cxxf77blas {

template <typename ENUM>
    typename std::enable_if<std::is_same<ENUM,Transpose>::value, char>::type
    getF77BlasChar(ENUM trans);

template <typename ENUM>
    typename std::enable_if<std::is_same<ENUM,Diag>::value, char>::type
    getF77BlasChar(ENUM diag);

template <typename ENUM>
    typename std::enable_if<std::is_same<ENUM,StorageUpLo>::value, char>::type
    getF77BlasChar(ENUM upLo);

} } // namespace cxxf77blas, flens

#endif // FLENS_AUXILIARY_CXXF77BLAS_H
