#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    GeMatrix<FullStorage<double> >   A(4,4);
    DenseVector<Array<double> >      x(4), y(4);

    A =  1,  2,  3,  4,
         5,  6,  7,  8,
         9, 10, 11, 12,
        13, 14, 15, 16;

    x =  1,  2,  3,  4;

    cout << "A = " << A << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    int m = A.numRows();
    int n = A.numCols();

//
//  compute  y = A*x
//
    cxxblas::gemv(m, n, 1.0, false, false,
                  A.data(), A.strideRow(), A.strideCol(),
                  x.data(), x.stride(),
                  0.0,
                  y.data(), y.stride());

    cout << "y = A*x = " << y << endl;

//
//  compute  y = A^T*x
//
    cxxblas::gemv(m, n, 1.0, true, false,
                  A.data(), A.strideRow(), A.strideCol(),
                  x.data(), x.stride(),
                  0.0,
                  y.data(), y.stride());

    cout << "y = A^T*x = " << y << endl;

//
//  compute  y = 1.5*y + A^T*x
//
    cxxblas::gemv(m, n, 1.0, true, false,
                  A.data(), A.strideRow(), A.strideCol(),
                  x.data(), x.stride(),
                  1.5,
                  y.data(), y.stride());

    cout << "y = 1.5*y + A^T*x = " << y << endl;
}
