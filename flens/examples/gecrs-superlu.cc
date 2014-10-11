#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

///
/// Include the SuperLu stuff for double precision.
///
#include <slu_ddefs.h>

///
/// Again we hack a quick and dirty wrapper for SuperLU.  This time we call the
/// SuperLU expert driver routine *dgssvx*: As our sparse matrix comes in
/// compressed row storage format but SuperLU expects the compressed col
/// storage format we actually habe to solve $A^T X = B$ instead of $AX = B$.
/// This operation is not available when using $dgssvx'.  The wrapper is based
/// on `dlinsolx.c` from the SuperLU examples directory.
///

template <typename MA, typename PR, typename PC, typename VB>
int
dgssvx(GeCRSMatrix<MA>  &A,
       DenseVector<PC>  &pc,
       DenseVector<PR>  &pr,
       DenseVector<VB>  &b)
{
///
/// We keep it as simple as possible and assume $A$ is square.
///
    ASSERT(A.numRows()==A.numCols());
    ASSERT(b.length()==A.numRows());
    const int n = A.numRows();

///
/// Some *temporary* vectors for intermediate results.
///
    typedef typename DenseVector<PC>::NoView  TempIntVec;
    typedef typename DenseVector<VB>::NoView  TempDoubleVec;
    TempDoubleVec   x(n), r(n), c(n), ferr(1), berr(1);
    TempIntVec      etree(n);

    superlu_options_t   options;
    SuperLUStat_t       stat;
    SuperMatrix         _A, _L, _U, _B, _X;

///
/// Set the transpose flag so that we solve $A^T X = B$.  Check the SuperLU
/// documentation for further details.
///
    set_default_options(&options);
    options.Equil           = YES;
    options.Trans           = TRANS;
    options.DiagPivotThresh = double(1);
    // Add more functionalities than the defaults.
    options.PivotGrowth = YES;        // Compute reciprocal pivot growth
    options.ConditionNumber = YES;    // Compute reciprocal condition number
    options.IterRefine = SLU_DOUBLE;  // Perform double-precision refinement

///
/// Because `A` has *CRS* format first pass `A.engine().cols().data()`and then
/// `A.engine().rows().data()`.  After creation `A_` holds a reference of $A^T$.
///
    dCreate_CompCol_Matrix(&_A, n, n, A.engine().numNonZeros(),
                           A.engine().values().data(),
                           A.engine().cols().data(),
                           A.engine().rows().data(),
                           SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&_B, n, 1, b.data(), b.length(),
                         SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&_X,
                         x.length(), 1, x.data(), x.length(),
                         SLU_DN, SLU_D, SLU_GE);

    StatInit(&stat);
    int             info;
    char            equed;
    int             lwork = 0; // let the solver allocate its own memory
    void            *work = 0;
    double          rpg, rcond;
    mem_usage_t     mem_usage;

///
/// Solve the system and compute the condition number and error bounds using
/// dgssvx.  Except for the computed solution this trashy SuperLU wrapper
/// ignores most of these computations.  However, in a real application some
/// of this stuff could be really useful (e.g. `rcond`!).
///
    dgssvx(&options, &_A, pc.data(), pr.data(), etree.data(), &equed,
           r.data(), c.data(), &_L, &_U, work, lwork, &_B, &_X,
           &rpg, &rcond, ferr.data(), berr.data(), &mem_usage, &stat, &info);

///
/// Overwrite $b$ with the solution $x$. (Quick and dirty, remember?)
///
    b = x;

    Destroy_SuperMatrix_Store(&_A);
    Destroy_SuperMatrix_Store(&_B);
    Destroy_SuperNode_Matrix(&_L);
    Destroy_CompCol_Matrix(&_U);
    StatFree(&stat);
    return info;
}

int
main()
{
///
/// SuperLU requires an index base of zero.  Here we set the default index base
/// of the storage scheme to zero via `IndexBaseZero`.  Note that the template
/// *Coord* contains the *CoordRowColCmp* class as third parameter.  That's
/// because we want to convert the storage to *compressed row storage* later.
///
    typedef int                                              IndexType;
    typedef IndexBaseZero<IndexType>                         IndexBase;
    typedef CoordStorage<double, CoordRowColCmp, IndexBase>  Coord;

    const IndexType m = 5;
    const IndexType n = 5;
    GeCoordMatrix<Coord>  A_(m, n);

///
/// Again we setup the matrix from the SuperLU user guide.  So here the values.
///
    const double s = 19,
                 u = 21,
                 p = 16,
                 e =  5,
                 r = 18,
                 l = 12;
    A_(0,0) += s;
    A_(1,1) += u;
    A_(2,2) += p;
    A_(3,3) += e;
    A_(4,4) += r;
    A_(1,0) += l; A_(2,1) += l; A_(4,0) += l; A_(4,1) += l;
    A_(0,2) += u; A_(0,3) += u; A_(3,4) += u;


///
/// Convert to *compressed row storage*.
///
    GeCRSMatrix<CRS<double, IndexBase> >  A = A_;


///
/// Just for curiosity: Compare the two formats.
///
    std::cout << "A_ = " << A_ << endl;
    std::cout << "A = " << A << endl;

///
/// Setup the right-hand side $b$.  We set all values to $1$.
///
    DenseVector<Array<double, IndexBase> >   b(m);
    b = 1;

///
/// Call the SuperLU solver for solving $Ax = b$.  Note that on exit $b$ is
/// overwritten with the solution $x$.
///
    DenseVector<Array<IndexType, IndexBase> >  pr(m), pc(n);
    int info = dgssvx(A, pr, pc, b);

///
/// Check the info code
///
    if (info==0) {
        cout << "x = " << b << endl;
    } else {
        cout << "SuperLU dgssv:  info = " << info << endl;
    }

    return info;
}
