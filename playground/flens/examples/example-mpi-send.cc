#include <cxxstd/iostream.h>

#define USE_PLAYGROUND
#define WITH_MPI
#include <flens/flens.cxx>

using namespace std;
using namespace flens;
using namespace mpi;

typedef double   T;

int
main(int argc, char* argv[])
{
    ///
    ///  Define some convenient typedefs for the vector types
    ///
    typedef DenseVector<Array<T> >              Vector;
    typedef Vector::IndexType                   IndexType;

    const IndexType m = 3;
    Vector x(m);

    ///
    /// Inititialize MPI enviroment
    ///
    MPI_init(argc, argv);

    ///
    /// Get rank
    ///
    int rank = MPI_rank();

    if ( rank==0 ) {
        ///
        /// Receive vectors and print them
        ///
        for (int i=1; i<MPI_size(); ++i) {
            MPI_recv(x, i);
            cout << "Rank " << i << " sent: ";
            cout << "x = " << x << endl;
        }

    } else {
        ///
        /// fill vector with numbers
        ///
        x = T(-rank), T(0), T(rank);

        ///
        /// Send vector to rank 0
        /// works also with GeMatrix
        ///
        MPI_send(x, 0);
    }
    ///
    /// Finalize MPI environment
    ///
    MPI_finalize();

    return 0;
}
