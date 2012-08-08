///
/// Define `USE_CXXLAPACK` before you include the FLENS headers
///
#ifndef USE_CXXLAPACK
#define USE_CXXLAPACK
#endif
#include <flens/flens.cxx>

using namespace flens;

///
/// Our high-level interface gets its own namespace
///
namespace mylapack {

///
/// We define our high-level interface for `getrf` it simply calls the 
/// CXXLAPACK interface `getrf` for the LAPACK functions `dgetrf`/`zgetrf`.
///
template <typename MA, typename VPIV>
typename GeMatrix<MA>::IndexType
trf(GeMatrix<MA> &A, DenseVector<VPIV> &piv)
{
    typedef typename GeMatrix<MA>::IndexType IndexType;
    return cxxlapack::getrf<IndexType>(A.numRows(),
                                       A.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       piv.data());
}

///
/// Analogously we define a high-level interface for `getrs`.  Note that
/// here the right hand side `b` is a vector.  So you might want to implement
/// another variant of this function were the right hand side is a general
/// matrix
///
template <typename MA, typename VPIV, typename VB>
void
trs(Transpose trans, const GeMatrix<MA> &A, const DenseVector<VPIV> &piv,
    DenseVector<VB> &b)
{
    typedef typename GeMatrix<MA>::IndexType IndexType;
    cxxlapack::getrs<IndexType>(lapack::getF77Char(trans),
                                A.numRows(),
                                1,
                                A.data(),
                                A.leadingDimension(),
                                piv.data(),
                                b.data(),
                                b.length());
}

} // namespace mylapack
