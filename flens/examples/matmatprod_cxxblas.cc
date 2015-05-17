#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    GeMatrix<FullStorage<double> >   A(4,4), B(4,4), C(4,4);

    A =  1,  2,  3,  4,
         5,  6,  7,  8,
         9, 10, 11, 12,
        13, 14, 15, 16;

    B = 17, 18, 19, 20,
        21, 22, 23, 24,
        25, 26, 27, 28,
        29, 30, 31, 32;

    auto U = A.upper();
    auto S = A.upper().symmetric();

    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "S = " << S << endl;
    cout << "U = " << U << endl;

//
//  compute  C = A*B
//
    cxxblas::gemm(C.numRows(), C.numCols(), A.numCols(),
                  1.0,
                  false, false, A.data(), A.strideRow(), A.strideCol(),
                  false, false, B.data(), B.strideRow(), B.strideCol(),
                  0.0,
                  C.data(), C.strideRow(), C.strideCol());

    cout << "C = A*B = " << C << endl;

//
//  compute  C = A^T*B
//
    cxxblas::gemm(C.numRows(), C.numCols(), A.numCols(),
                  1.0,
                  true,  false, A.data(), A.strideRow(), A.strideCol(),
                  false, false, B.data(), B.strideRow(), B.strideCol(),
                  0.0,
                  C.data(), C.strideRow(), C.strideCol());

    cout << "C = A^T*B = " << C << endl;

//
//  compute  C = S*B
//
    cxxblas::symm(true, C.numRows(), C.numCols(),
                  1.0,
                  (S.upLo()==Lower),
                  S.data(), S.strideRow(), S.strideCol(),
                  B.data(), B.strideRow(), B.strideCol(),
                  0.0,
                  C.data(), C.strideRow(), C.strideCol());

    cout << "C = S*B = " << C << endl;

//
//  compute  B = U*B
//
    cxxblas::trmm(true, C.numRows(), C.numCols(),
                  1.0,
                  (U.upLo()==Lower), false, false, (U.diag()==Unit),
                  U.data(), U.strideRow(), U.strideCol(),
                  B.data(), B.strideRow(), B.strideCol());

    cout << "B = U*B = " << B << endl;
}
