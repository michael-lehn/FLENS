#include <iostream>
///
///  With header __flens.cxx__ all of FLENS gets included.
///
///  :links:  __flens.cxx__ -> file:flens/flens.cxx
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

typedef double   T;

int
main()
{
    ///
    ///  Define some convenient typedefs for the matrix/vector types
    ///  of our system of linear equations.
    ///
    typedef GeMatrix<FullStorage<T, ColMajor> >   Matrix;
    typedef DenseVector<Array<T> >                Vector;

    const int n = 5;

    Matrix   A(n, n), VL(n, n), VR(n, n);
    Vector   wr(n), wi(n);


    A =  2,   3,  -1,   0,  2,
        -6,  -5,   0,   2, -6,
         2,  -5,   6,  -6,  2,
         2,   3,  -1,   0,  8,
        -6,  -5,  10,   2, -6;

    cerr << "A = " << A << endl;

    ///
    ///  Vector for workspace.  If this vector has zero length then
    ///  __lapack::ev__ will do a worksize query and also resize `work`.
    ///
    Vector   work;

    ///
    ///  You also could do a worksize query manually
    ///
    // int     optSize = ev_wsq(true, true, A);
    // Vector  work(optSize);

    ///
    ///  Call __lapack::ev__ to compute eigenvalues $w = w_r+i w_i$,
    ///  left eigenvectors $V_L$ and right eigenvectors $V_R$.
    ///
    ///  :links: __lapack::ev__ -> file:flens/lapack/eig/ev.h
    lapack::ev(true, true, A, wr, wi, VL, VR, work);

    cerr << "wr = " << wr << endl;
    cerr << "wi = " << wi << endl;
    cerr << "VL = " << VL << endl;
    cerr << "VR = " << VR << endl;
}
