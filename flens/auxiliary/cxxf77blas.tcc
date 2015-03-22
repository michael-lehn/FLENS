#ifndef FLENS_AUXILIARY_CXXF77BLAS_TCC
#define FLENS_AUXILIARY_CXXF77BLAS_TCC 1

namespace flens { namespace cxxf77blas {

template <typename ENUM>
typename std::enable_if<std::is_same<ENUM,Transpose>::value, char>::type
getF77BlasChar(ENUM trans)
{
    if (trans==NoTrans) {
        return 'N';
    } else if (trans==Trans) {
        return 'T';
    } else if (trans==Conj) {
        return 'R';
    } else if (trans==ConjTrans) {
        return 'C';
    } else {
        ASSERT(0);
        return '?';
    }
}

template <typename ENUM>
typename std::enable_if<std::is_same<ENUM,Diag>::value, char>::type
getF77BlasChar(ENUM diag)
{
    if (diag==Unit) {
        return 'U';
    }
    return 'N';
}

template <typename ENUM>
typename std::enable_if<std::is_same<ENUM,StorageUpLo>::value, char>::type
getF77BlasChar(ENUM upLo)
{
    if (upLo==Upper) {
        return 'U';
    }
    return 'L';
}

} } // namespace cxxf77blas, flens

#endif // FLENS_AUXILIARY_CXXF77BLAS_TCC
