#include <iostream>
///
///  With header __flens.cxx__ all of FLENS gets included.
///
///  :links:  __flens.cxx__ -> file:flens/flens.cxx
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

typedef double   T;

namespace flens { namespace mylapack {

template <typename MA, typename VPIV>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
trf(MA &&A, VPIV &&piv)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename RemoveRef<VPIV>::Type  VectorPiv;

//
//  Create views of the arguments
//
    typename MatrixA::View    _A      = A;
    typename VectorPiv::View  _piv    = piv;

//
//  Make the views one-based
//
    _A.changeIndexBase(1,1);
    _piv.changeIndexBase(1);

    IndexType info = lapack::trf(_A, _piv);

    const IndexType diff = piv.firstIndex() - _piv.firstIndex();

    for (IndexType i=1; i<=_piv.length(); ++i) {
        _piv(i) += diff;
    }

    return info;
}

} } // namespace mylapack, flens

int
main()
{
    ///
    ///  Define some convenient typedefs for the matrix type
    ///  of our system of linear equations.
    ///
    typedef GeMatrix<FullStorage<T> >           Matrix;

    ///
    ///  We also need an extra vector type for the pivots.  The type of the
    ///  pivots is taken for the system matrix.
    ///
    typedef Matrix::IndexType                   IndexType;
    typedef DenseVector<Array<IndexType> >      IndexVector;

    ///
    ///  Set up the baby problem ...
    ///
    const IndexType m = 4,
                    n = 5;

    ///
    /// Zero based matrix and vector
    ///
    Matrix            Ab(m, n, 0, 0);
    IndexVector       piv(m, 0);

    Ab =  2,   3,  -1,   0,  20,
         -6,  -5,   0,   2, -33,
          2,  -5,   6,  -6, -43,
          4,   6,   2,  -3,  49;

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
