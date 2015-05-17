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
}

template <typename MA, typename VV>
void
eigenvalues(const GeMatrix<MA> &A, DenseVector<VV> &v)
{
    DenseVector<VV> ws;

    eigenvalues(A, v, ws);
}

int
main()
{
    typedef GeMatrix<FullStorage<double> >  DGeMatrix;
    typedef DenseVector<Array<double> >     DDenseVector;

    // Use case 1: Use optimal workspace (performance for memory).  We have
    //             a typical case where the same problem size is involved.  So
    //             the workspace can be reused.
    {
        DGeMatrix     A(100, 100);
        DDenseVector  v(100);

        auto workspace = eigenvalues_wsq(A);
        DDenseVector  ws(workspace.second);

        for (int i=0; i<666; ++i) {
            fillRandom(A);
            eigenvalues(A, v, ws);
        }
    }

    // Use case 2: Use minimal workspace (memory for performance)
    {
        DGeMatrix     A(100, 100);
        DDenseVector  v(100);

        auto workspace = eigenvalues_wsq(A);
        DDenseVector  ws(workspace.first);

        for (int i=0; i<666; ++i) {
            fillRandom(A);
            eigenvalues(A, v, ws);
        }
    }

    // Use case 3: Allocate (optimal) workspace on first call.
    {
        DGeMatrix     A(100, 100);
        DDenseVector  v(100);
        DDenseVector  ws;

        for (int i=0; i<666; ++i) {
            fillRandom(A);
            eigenvalues(A, v, ws);
        }
    }

    // NO-Use case 4: Workspace gets allocated for each call!!!
    {
        DGeMatrix     A(100, 100);
        DDenseVector  v(100);

        for (int i=0; i<666; ++i) {
            fillRandom(A);
            eigenvalues(A, v);
        }
    }

    // Use case 5: Ok, workspace is needed only once.
    {
        DGeMatrix     A(100, 100);
        DDenseVector  v(100);

        fillRandom(A);
        eigenvalues(A, v);
    }
}
