#include <cassert>
#include <iostream>
#include <utility>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

template <typename MA>
std::pair<long, long>
eigenvalues_wsq(const GeMatrix<MA> &A)
{
    // Compute and return minimal and optimal workspace sizes
    return std::pair<long,long>(A.numRows()*4, A.numRows()*256);
}

template <typename MA, typename VV, typename VWS>
void
eigenvalues(const GeMatrix<MA> &A, DenseVector<VV> &v, DenseVector<VWS> &ws)
{
    auto wsq = eigenvalues_wsq(A);
    if (ws.length()==0) {
        ws.resize(wsq.second);
    }
    assert(ws.length()>=wsq.first);

    //
    // ... compute eigenvalues v of A ... here the actual work happens!!!
    //
    cout << "Computing eigenvalue "
         << "(Using a workspace of size " << ws.length() << ")"
         << endl;
}

template <typename MA, typename VV>
void
eigenvalues(const GeMatrix<MA> &A, DenseVector<VV> &v)
{
    DenseVector<VV> ws;

    eigenvalues(A, v, ws);
}


double buffer_[DIM*DIM + DIM + WORKSPACE];

int
main()
{
    typedef GeMatrix<FullStorage<double> >  DGeMatrix;
    typedef DenseVector<Array<double> >     DDenseVector;

    DGeMatrix::View    A  = DGeMatrix::EngineView(DIM, DIM,
                                                  &buffer_[0],
                                                  1, DIM);

    DDenseVector::View v  = DDenseVector::EngineView(DIM,
                                                     &buffer_[DIM*DIM]);

    DDenseVector::View ws = DDenseVector::EngineView(WORKSPACE,
                                                     &buffer_[DIM*DIM+DIM]);

#   if WORKSPACE==0
    // Perform workspace query and return

    auto workspace = eigenvalues_wsq(A);
    cout << workspace.first << " " << workspace.second << endl;

#   else
    // Solve the problem with given workspace

    for (int i=0; i<6; ++i) {
        fillRandom(A);
        eigenvalues(A, v, ws);
    }
#   endif
}
