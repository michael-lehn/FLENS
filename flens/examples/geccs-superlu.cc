#include <iostream>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

///
/// Include the SuperLu stuff for double precision.
///
#include <slu_ddefs.h>

///
/// Wrapper for the SuperLU solver dgssv.  Looks complicated but it's actually
/// pretty simple.  We create SuperLU matrix wrapper for our matrices and then
/// call the solver.  Anyway, in future we will hide the details in FLENS.
///
template <typename MA, typename PR, typename PC, typename VB>
int
dgssv(GeCCSMatrix<MA>  &A,
      DenseVector<PC>  &pc,
      DenseVector<PR>  &pr,
      DenseVector<VB>  &b)
{
    ASSERT(pr.length()==A.numRows());
    ASSERT(pc.length()==A.numCols());
    superlu_options_t   options;
    SuperLUStat_t       stat;
    SuperMatrix         _A, _L, _U, _B;
    dCreate_CompCol_Matrix(&_A,
                           A.numRows(), A.numCols(), A.engine().numNonZeros(),
                           A.engine().values().data(),
                           A.engine().rows().data(),
                           A.engine().cols().data(),
                           SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&_B,
                         b.length(), 1, b.data(), b.length(),
                         SLU_DN, SLU_D, SLU_GE);
    set_default_options(&options);
    options.ColPerm = NATURAL;
    StatInit(&stat);
    int info;
    dgssv(&options, &_A, pc.data(), pr.data(), &_L, &_U, &_B, &stat, &info);
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
/// of the storage scheme to zero via `IndexBaseZero`.
///
    typedef int                                              IndexType;
    typedef IndexBaseZero<IndexType>                         IndexBase;
    typedef CoordStorage<double, CoordColRowCmp, IndexBase>  Coord;

///
/// Alternative we could specify it through the constructor (see class API).
///
    const IndexType m = 5;
    const IndexType n = 5;
    GeCoordMatrix<Coord>  A_(m, n);

///
/// We setup the matrix from the SuperLU user guide.  So here the values.
///
    const double s = 19,
                 u = 21,
                 p = 16,
                 e =  5,
                 r = 18,
                 l = 12;
///
/// Matrices in coordinate storage are usually used in FEM for assembling
/// the stiffness matrix.  Values will be accumulated.  Note that an assignment
/// like 'A(2,3) = ' is not allowed.
///
    A_(0,0) += s;
    A_(1,1) += u;
    A_(2,2) += p;
    A_(3,3) += e;
    A_(4,4) += r;
    A_(1,0) += l; A_(2,1) += l; A_(4,0) += l; A_(4,1) += l;
    A_(0,2) += u; A_(0,3) += u; A_(3,4) += u;


///
/// Convert to compressed column storage.
///
    GeCCSMatrix<CCS<double, IndexBase> >  A = A_;


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
    int info = dgssv(A, pr, pc, b);

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
