#include <iostream>
/*
#define CXXBLAS_DEBUG_OUT(X)                                                \
   flens::verbose::ClosureLog::append(false) << "   => BLAS: " << X << ";"; \
*/


#include <flens/debug/aux/aux.h>
#include <flens/flens.cxx>
#include <flens/debug/aux/aux.tcc>

using namespace std;
using namespace flens;


int
main()
{
    double x = lapack::hex<double>("0xe2 0x92 0x46 0x92 0x1e 0x29 0x5d 0x83");
    double y = lapack::hex<double>("0x88 0xbb 0x88 0x11 0x5e 0x2d 0x2e 0x03");

    std::cerr.precision(1550);
    cerr << "x = " << x << endl;
    cerr << "y = " << y << endl;

    double result = lapack::lapy2(x, y);

    cerr << "result = " << result << endl;

    cerr << "hex(x) =      " << lapack::hex(x) << endl;
    cerr << "hex(y) =      " << lapack::hex(y) << endl;
    cerr << "hex(result) = " << lapack::hex(result) << endl;

}

