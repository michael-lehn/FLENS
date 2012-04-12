#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

#ifdef FLENS_DEBUG_CLOSURES
#   define CLOSURELOG_START(x) verbose::ClosureLog::start(x)
#   define CLOSURELOG_STOP     verbose::ClosureLog::stop()
#else
#   define CLOSURELOG_START(x)
#   define CLOSURELOG_STOP
#endif

int
main()
{
    GEMatrix A(3,3), B(3,3), C;
    A = 1, 2, 3,
        5, 6, 7,
        5, 4, 3;

    B = 9, 4, 2,
        5, 1, 9,
        6, 7, 1;

    CLOSURELOG_START("mylogfile");

    ///
    /// Compute $C = 2 A B$
    ///
    C = 2*A*B;

    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "C = " << C << endl;

    CLOSURELOG_STOP;

    return 0;
}

