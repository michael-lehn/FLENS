///
/// Define `USE_CXXLAPACK` before you include the FLENS headers
///
#define USE_CXXLAPACK
#include <flens/flens.cxx>

using namespace flens;

namespace mylapack {

///
/// The trait is of form `RestrictTo<Cond,B>::Type`.  If `Cond` is
/// true then `RestrictTo<Cond,B>::Type` equals `B`.  This will be the return
/// type.  If `Cond` is false then the trait can not be instantiated.  In this
/// case the function can not be instantiated either and the whole thing gets
/// ignored.
///
/// Note that we also need the `RemoveRef` trait (yeah, C++ brain fuck!).  If
/// `MA` is of type `T &` then `RemoveRef<T>::Type` is of type `T`.
///
template <typename MA, typename VPIV>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
trf(MA &&A, VPIV &&piv)
{
    typedef typename RemoveRef<MA>::Type::IndexType  IndexType;
    return cxxlapack::getrf<IndexType>(A.numRows(),
                                       A.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       piv.data());
}

///
/// For `trs` things are like above.  But now the return type is simpler and
/// for traits-newbies it might be easier to read.
///
template <typename MA, typename VPIV, typename VB>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsDenseVector<VB>::value,
         void>::Type
trs(Transpose trans, const MA &A, const VPIV &piv, VB &&b)
{
    typedef typename MA::IndexType IndexType;
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
