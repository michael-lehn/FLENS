using namespace flens;

namespace mylapack {

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
/// Note the `&&` in front of b.  This means that `b` can now be non-const
/// rvalue objects.  Such objects arise when we create views using the
/// underscore operator, e.g. by `b(_(1,n))`.
///
template <typename MA, typename VPIV, typename VB>
void
trs(Transpose trans, const GeMatrix<MA> &A, const DenseVector<VPIV> &piv,
    DenseVector<VB> &&b)
///
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
