#include <iostream>

///
/// Define `USE_CXXLAPACK` before you include the FLENS headers
///
#define USE_CXXLAPACK
#include <flens/flens.cxx>

using namespace std;
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
    return cxxlapack::getrf(A.numRows(),
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


typedef double   T;

int
main()
{
    typedef GeMatrix<FullStorage<T> >           Matrix;
    typedef DenseVector<Array<T> >              Vector;

    typedef Matrix::IndexType                   IndexType;
    typedef DenseVector<Array<IndexType> >      IndexVector;

    const IndexType n = 4;

    Matrix         A(n,n);
    Vector         b(n);
    IndexVector    piv(n);

    A =  2,   3,  -1,   0,
        -6,  -5,   0,   2,
         2,  -5,   6,  -6,
         4,   6,   2,  -3.1;

    b = 20,
       -33,
       -43,
        49;

    cerr << "A = " << A << endl;
    cerr << "b = " << b << endl;

///
/// We now can use our interface just the way we are using __FLENS-LAPACK__ ...
///
/// :links: __FLENS-LAPACK__ -> doc:flens/lapack/lapack
///
    IndexType info = mylapack::trf(A, piv);

    if (info==0) {
///
///     ... for a user only the namspace differs.
///
        mylapack::trs(NoTrans, A, piv, b);

        cerr << "x = " << b << endl;
    } else {
        cerr << "A is numerically singular." << endl;
    }
}

