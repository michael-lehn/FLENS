#include <iostream>
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
    Vector x(m), sum(m), max(m), min(m);

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
        /// fill vector with numbers
        ///
        x = T(1), T(2), T(3);

    } 
    
    ///
    /// Broadcast vector from rank 0, 
    /// works also with GeMatrix
    ///
    MPI_bcast(x, 0);
    
    ///
    /// Print result
    ///
    cout << "Rank = " << rank << ", x = " << x << endl;
    
    MPI_finalize();

    return 0;
}
