#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    typedef double                               T;
    typedef DenseVector<Array<T> >               DEVector;
    typedef GeMatrix<FullStorage<T, ColMajor> >  GEMatrix;

    DEVector x(3);
    x = 1, 2, 3;

    GEMatrix A(3,3);
    A = 1, 2, 3,
        5, 6, 7,
        5, 4, 3;

    DEVector y(3);
    ///
    /// We turn on logging of the closure evaluation, i.e. the evaluation
    /// of linear algebra expressions that are coded through overloaded
    /// operators.
    ///
    verbose::ClosureLog::start("mylogfile");

    ///
    /// Compute $y = A x$
    ///
    y = -A*x;

    cout << "x = " << x << endl;
    cout << "A = " << A << endl;
    cout << "y = " << y << endl;

    //y = conjugate(transpose(A))*x;
    y = conjTrans(A)*x;
    cout << "x = " << x << endl;
    cout << "A = " << A << endl;
    cout << "y = " << y << endl;

    ///
    /// Stop logging.
    ///
    verbose::ClosureLog::stop();

    return 0;
}
