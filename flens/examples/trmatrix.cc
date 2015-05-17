#include <flens/flens.cxx>
#include <algorithm>
#include <iostream>

using namespace std;
using namespace flens;


// Initialize a triangular/trapezoidal matrix column wise with 1, 2, 3, ...
template <typename MA>
void
myInit(TrMatrix<MA> &A)
{
    int counter = 1;

    int m = A.numRows();
    int n = A.numCols();

    if (A.upLo()==Upper) {
        if (A.diag()==NonUnit) {
            for (int j=1; j<=n; ++j) {
                for (int i=1; i<=std::min(j,m); ++i) {
                    A(i,j) = counter++;
                }
            }
        } else {
            for (int j=1; j<=n; ++j) {
                for (int i=1; i<std::min(j,m); ++i) {
                    A(i,j) = counter++;
                }
            }
        }
    } else {
        if (A.diag()==NonUnit) {
            for (int j=1; j<=n; ++j) {
                for (int i=j; i<=m; ++i) {
                    A(i,j) = counter++;
                }
            }
        } else {
            for (int j=1; j<=n; ++j) {
                for (int i=j+1; i<=m; ++i) {
                    A(i,j) = counter++;
                }
            }
        }
    }
}

int
main()
{
    TrMatrix<FullStorage<double> >  A(3, Upper);
    TrMatrix<FullStorage<double> >  B(3, Upper, Unit);
    TrMatrix<FullStorage<double> >  C(3, 4, Upper);
    TrMatrix<FullStorage<double> >  D(3, 4, Upper, Unit);
    TrMatrix<FullStorage<double> >  E(4, 3, Upper);
    TrMatrix<FullStorage<double> >  F(4, 3, Upper, Unit);

    TrMatrix<FullStorage<double> >  G(3, Lower);
    TrMatrix<FullStorage<double> >  H(3, Lower, Unit);
    TrMatrix<FullStorage<double> >  I(3, 4, Lower);
    TrMatrix<FullStorage<double> >  J(3, 4, Lower, Unit);
    TrMatrix<FullStorage<double> >  K(4, 3, Lower);
    TrMatrix<FullStorage<double> >  L(4, 3, Lower, Unit);

    myInit(A);
    myInit(B);
    myInit(C);
    myInit(D);
    myInit(E);
    myInit(F);

    myInit(G);
    myInit(H);
    myInit(I);
    myInit(J);
    myInit(K);
    myInit(L);

    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "C = " << C << endl;
    cout << "D = " << D << endl;
    cout << "E = " << E << endl;
    cout << "F = " << F << endl;

    cout << "G = " << G << endl;
    cout << "H = " << H << endl;
    cout << "I = " << I << endl;
    cout << "J = " << J << endl;
    cout << "K = " << K << endl;
    cout << "L = " << L << endl;
}
