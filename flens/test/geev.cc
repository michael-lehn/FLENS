//
// compile with:
//
// clang++ -std=c++0x -Wall -I../../ geev.cc
//
//

#include <iostream>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

// typedef mpf_class   T;
typedef long double   T;

int
main()
{
    typedef GeMatrix<FullStorage<T> >   GeMatrix;
    typedef GeMatrix::IndexType         IndexType;
    typedef DenseVector<Array<T> >      DenseVector;

    const Underscore<IndexType> _;

    const IndexType n = 3;

    GeMatrix            A(n, n), VL(n,n), VR(n,n);
    DenseVector         wr(n), wi(n), work;

    A = 1, 2, 3,
        4, 5, 6,
        7, 8, 9;

    cerr << "A = " << A << endl;

    lapack::ev(true, true, A, wr, wi, VL, VR, work);


    cerr << "-> A = " << A << endl;
    cerr << "-> wr = " << wr << endl;
    cerr << "-> wi = " << wi << endl;
    cerr << "-> VL = " << VL << endl;
    cerr << "-> VR = " << VR << endl;
}
