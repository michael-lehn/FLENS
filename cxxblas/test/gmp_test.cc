// clang++ gmp_test.cc -I../../ -I/opt/local/include/ -L /opt/local/lib -lgmpxx -lgmp


#include <gmpxx.h>
#include <iostream>
#include <cxxblas/cxxblas.cxx>

using namespace std;

typedef mpq_class   T;

int
main()
{
    const int n = 30;
    T a[n], b[n], c;
    
    for (int i=0; i<n; ++i) {
        a[i] = T(n+i)/n;
        
        b[i] = 1/T(n-i);
    }

    cxxblas::dot(n, a, 1, b, 1, c);

    cout << "a*b = " << c << "\n";
    
    cxxblas::axpy(n, c, b, 1, a, 1);
    for (int i=0; i<n; ++i) {
        cout << i << ") " << a[i] << endl;
    }
}