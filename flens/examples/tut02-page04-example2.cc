#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef double                               T;
    typedef DenseVector<Array<T> >               DEVector;
    typedef GeMatrix<FullStorage<T, ColMajor> >  GEMatrix;

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
    /// Compute $v = x + (y + z)$ note that this *must not* be evaluated
    /// as $v = (x + y) + z
    ///
    DEVector v = x + (y + z);

    ///
    /// Compute $z = A (x + y)$ note that this *must not* be evaluated
    /// as $z = Ax + Ay$.
    ///
    z = A*(x+y);

    cout << "x = " << x << endl;

    ///
    /// Stop logging.
    ///
    verbose::ClosureLog::stop();

    return 0;
}
