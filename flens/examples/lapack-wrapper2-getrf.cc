#include <cxxstd/iostream.h>
///
///  With header __flens.cxx__ all of FLENS gets included.
///
///  :links:  __flens.cxx__ -> file:flens/flens.cxx
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

typedef double   T;

//
//  Traits for creating one based matrix/vector views
//
namespace flens {

template <typename A>
struct OneBased
{
};

template <typename MA>
struct OneBased<GeMatrix<MA> >
{
    typedef typename MA::ElementType      T;
    typedef typename MA::IndexType        IndexType;
    typedef IndexOptions<IndexType, 1>    IndexBase;

    static const StorageOrder order = MA::order;

    typedef GeMatrix<FullStorageView<T, order, IndexBase> >         View;
    typedef GeMatrix<ConstFullStorageView<T, order, IndexBase> >    ConstView;
};

template <typename VX>
struct OneBased<DenseVector<VX> >
{
    typedef typename VX::ElementType      T;
    typedef typename VX::IndexType        IndexType;
    typedef IndexOptions<IndexType, 1>    IndexBase;

    typedef DenseVector<ArrayView<T, IndexBase> >         View;
    typedef DenseVector<ConstArrayView<T, IndexBase> >    ConstView;
};

} // namespace flens

//
//  LAPACK wrapper for non-one based indices
//
namespace flens { namespace mylapack {

template <typename MA, typename VPIV>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
trf(MA &&A, VPIV &&piv)
{
///
/// Remove references from rvalue types
///
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename RemoveRef<VPIV>::Type  VectorPiv;

///
/// Create views of the arguments
///
    typename OneBased<MatrixA>::View    A_    = A;
    typename OneBased<VectorPiv>::View  piv_  = piv;

    A_.changeIndexBase(1,1);
    piv_.changeIndexBase(1);

///
/// Make the views one-based
///
    IndexType info = lapack::trf(A_, piv_);

    const IndexType diff = piv.firstIndex() - piv_.firstIndex();

    for (IndexType i=1; i<=piv_.length(); ++i) {
        piv_(i) += diff;
    }

    return info;
}

} } // namespace mylapack, flens

int
main()
{
    ///
    ///  Define an index-option type for zero based indices
    ///
    typedef IndexOptions<int, 0>                           ZeroBased;

    ///
    ///  Define some convenient typedefs for the matrix type
    ///  of our system of linear equations.
    ///
    typedef GeMatrix<FullStorage<T, ColMajor, ZeroBased> > Matrix;

    ///
    ///  We also need an extra vector type for the pivots.  The type of the
    ///  pivots is taken for the system matrix.
    ///
    typedef Matrix::IndexType                              IndexType;
    typedef DenseVector<Array<IndexType, ZeroBased> >      IndexVector;

    ///
    ///  Set up the baby problem ...
    ///
    const IndexType m = 4,
                    n = 4;

    ///
    /// Zero based matrix and vector
    ///
    Matrix            Ab(m, n);
    IndexVector       piv(m);

    Ab = 2,  4,  4,  4,
         4,  2,  4,  4,
         4,  4,  2,  4,
         4,  4,  4,  2;

    cout << "Ab = " << Ab << endl;

    ///
    /// Compute the $LU$ factorization with __lapack::trf__
    ///
    mylapack::trf(Ab, piv);

    cout << "Ab.firstRow() = " << Ab.firstRow() << endl;
    cout << "Ab.firstCol() = " << Ab.firstCol() << endl;

    cout << "Ab = " << Ab << endl;
    cout << "piv = " << piv << endl;
}
