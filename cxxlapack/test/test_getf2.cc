// clang++ test_getf2.cc -I../../ -I/opt/local/include/ -L /opt/local/lib -lgmpxx -lgmp


#include <gmpxx.h>
#include <iostream>
#include <cxxblas/cxxblas.cxx>
#include <cxxlapack/cxxlapack.cxx>

using namespace std;

typedef mpq_class   T;

// typedef double   T;

int
main()
{
    const int m = 4;
    const int n = 5;
    T A[n][m];
    const int ldA = m;

    int iPiv[m];

    /*
    A[0][0]=  2;  A[1][0]=  3;   A[2][0]= -1;   A[3][0]=  0;   A[4][0]=  20;
    A[0][1]= -6;  A[1][1]= -5;   A[2][1]=  0;   A[3][1]=  2;   A[4][1]= -33;
    A[0][2]=  2;  A[1][2]= -5;   A[2][2]=  6;   A[3][2]= -6;   A[4][2]= -43;
    A[0][3]=  4;  A[1][3]=  6;   A[2][3]=  2;   A[3][3]= -3;   A[4][3]=  49;
    */

    A[0][0]=  3;  A[1][0]= -1;   A[2][0]=  0;   A[3][0]=  0;   A[4][0]=  4;
    A[0][1]= -1;  A[1][1]=  3;   A[2][1]= -1;   A[3][1]=  0;   A[4][1]= -2;
    A[0][2]=  0;  A[1][2]= -1;   A[2][2]=  3;   A[3][2]= -1;   A[4][2]= -7;
    A[0][3]=  0;  A[1][3]=  0;   A[2][3]= -1;   A[3][3]=  3;   A[4][3]=  8;

    for (int i=0; i<m; ++i) {
        for (int j=0; j<n; ++j) {
            cout.width(6);
            cout << A[j][i] << "  ";
        }
        cout << endl;
    }
    cout << endl;

    T *_A = &(A[0][0]);

    cxxlapack::gesv(cxxblas::ColMajor, m, n-m, _A, ldA, iPiv, _A + m*ldA, ldA);

    for (int i=0; i<m; ++i) {
        for (int j=0; j<n; ++j) {
            cout.width(6);
            cout << A[j][i] << "  ";
        }
        cout << endl;
    }
    cout << endl;
    for (int i=0; i<m; ++i) {
        cout.width(6);
        cout << iPiv[i] << "  ";
    }
    cout << endl << endl;
}
