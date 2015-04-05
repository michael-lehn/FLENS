#ifndef FLENS_AUXILIARY_CXXF77BLAS_TCC
#define FLENS_AUXILIARY_CXXF77BLAS_TCC 1

#include <flens/auxiliary/cxxf77blas.h>
#include <flens/auxiliary/macros.h>

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

//------------------------------------------------------------------------------
template <typename ENUM>
typename std::enable_if<std::is_same<ENUM,Transpose>::value,
                        Transpose>::type
getCxxBlasEnum(char trans)
{
    if ((trans=='N') || (trans=='n')) {
        return NoTrans;
    } else if ((trans=='T') || (trans=='t')) {
        return Trans;
    } else if ((trans=='C') || (trans=='c')) {
        return ConjTrans;
    } else if ((trans=='R') || (trans=='r')) {
        return Conj;
    }
    ASSERT(0);
    return NoTrans;
}

template <typename ENUM>
typename std::enable_if<std::is_same<ENUM,Diag>::value,
                        Diag>::type
getCxxBlasEnum(char diag)
{
    if (diag=='U') {
        return Unit;
    }
    ASSERT(diag=='N');
    return NonUnit;
}

template <typename ENUM>
typename std::enable_if<std::is_same<ENUM,StorageUpLo>::value,
                        StorageUpLo>::type
getCxxBlasEnum(char upLo)
{
    if (upLo=='U') {
        return Upper;
    }
    ASSERT(upLo=='L');
    return Lower;
}

} } // namespace cxxf77blas, flens

#endif // FLENS_AUXILIARY_CXXF77BLAS_TCC
