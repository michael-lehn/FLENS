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

    const T  alpha = 1.5,
             beta = 2.5;

    DEVector x(3), y(3), z(3);
    x = 1, 2, 3;
    y = 2, 3, 4;
    z = 3, 4, 5;

    GEMatrix A(3,3);
    A = 1, 2, 3,
        5, 6, 7,
        5, 4, 3;

    ///
    /// We turn on logging of the closure evaluation, i.e. the evaluation
    /// of linear algebra expressions that are coded through overloaded
    /// operators.
    ///
    verbose::ClosureLog::start("mylogfile");

    ///
    /// Compute $y = \beta y + \alpha A^T x + z$
    ///
    y = beta*y + alpha*transpose(A)*x + z;

    cout << "y = " << y << endl;

    ///
    /// Stop logging.
    ///
    verbose::ClosureLog::stop();

    return 0;
}
