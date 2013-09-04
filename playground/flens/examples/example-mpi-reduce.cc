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
    
    ///
    /// fill vector with numbers
    ///
    x = T(-rank), T(0), T(rank);
    
    ///
    /// Sum up all vectors x and store result at rank 0
    /// works also with GeMatrix
    ///
    MPI_reduce_sum(x, sum, 0);

    if ( rank==0 ) { 
        ///
        /// Rank 0 print result
        ///
        cout << "sum = " << sum << endl;

    } 
    
    MPI_finalize();

    return 0;
}
