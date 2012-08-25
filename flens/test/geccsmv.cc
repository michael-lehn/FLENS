#include <iostream>
#include <flens/flens.cxx>

#ifndef SEED
#define SEED  0
#endif

#ifndef MAX_M
#define MAX_M  500
#endif

#ifndef MAX_N
#define MAX_N  500
#endif

#ifndef MAX_NNZ
#define MAX_NNZ  3*MAX_M
#endif


using namespace flens;
using namespace std;

//
//  Create sparse matrix A and control matrix A_
//
template <typename CCS, typename FS>
void
setup(int m, int n, int max_nnz, int indexBase,
      GeCCSMatrix<CCS> &A, GeMatrix<FS> &A_)
{
    typedef typename GeCCSMatrix<CCS>::ElementType  ElementType;
    typedef CoordStorage<double, CoordColRowCmp>    Coord;

    const ElementType  Zero(0);

    A_.resize(m, n, indexBase, indexBase);
    A_ = Zero;

    //
    //  We first setup the sparse matrix B in coordinate storage.  Later we
    //  convert it to compressed row storage.
    //
    GeCoordMatrix<Coord>            B(m, n, 1, indexBase);


    for (int k=1; k<=max_nnz; ++k) {
        const int i = indexBase + rand() % m;
        const int j = indexBase + rand() % n;
        const int v1 = rand() % 10;
        const int v2 = rand() % 10;

        B(i,j) += v1;
        A_(i,j) += v1;

        B(i,j) -= v2;
        A_(i,j) -= v2;
    }

    //
    //  Convert coordinate storage matrix B to compressed row storage matrix A
    //
    A = B;

    //
    //  Convert compressed row storage matrix A to full storage matrix A__
    //
    typename GeMatrix<FS>::NoView  A__ = A;

    if (! lapack::isIdentical(A_, A__, "A_", "A__")) {
        cerr << "m =       " << m << endl;
        cerr << "n =       " << n << endl;
        cerr << "max_nnz = " << max_nnz << endl;
        ASSERT(0);
    }
}

template <typename CCS, typename FS>
void
mv(int m, int n, int max_nnz, const GeCCSMatrix<CCS> &A, const GeMatrix<FS> &A_)
{
    typedef typename GeCCSMatrix<CCS>::ElementType  ElementType;

    DenseVector<Array<ElementType> >  x(n), y, y_;

    for (int j=1; j<=n; ++j) {
        x(j) = rand() % 10;
    }

    y  = A  * x;
    y_ = A_ * x;

    if (! lapack::isIdentical(y, y_, "y", "y_")) {
        cerr << endl << "failed: y = A*x" << endl;
        ASSERT(0);
    }

    y  += A  * x;
    y_ += A_ * x;

    if (! lapack::isIdentical(y, y_, "y", "y_")) {
        cerr << endl << "failed: y += A*x" << endl;
        ASSERT(0);
    }

    y  -= A  * x;
    y_ -= A_ * x;

    if (! lapack::isIdentical(y, y_, "y", "y_")) {
        cerr << endl << "failed: y -= A*x" << endl;
        ASSERT(0);
    }

    ElementType  alpha, beta;
    for (int test=1; test<=20; ++test) {

        alpha = std::pow(2, 5 - std::max(1, rand() % 10));
        beta  = std::pow(2, 5 - std::max(1, rand() % 10));

        //
        // Reset y (and y_) to some random vector
        //
        for (int i=1; i<=m; ++i) {
            y(i) = rand() % 1000;
        }
        y_ = y;

        y  = beta*y + alpha*A  * x;
        y_ = beta*y_ + alpha*A_ * x;

        if (! lapack::isIdentical(y, y_, "y", "y_")) {
            cerr << endl << "failed: y = beta*y + alpha*A*x" << endl;
            cout << "alpha = " << alpha << endl;
            cout << "beta  = " << beta << endl;
            cout << "A_  = " << A_ << endl;
            cout << "x  = " << x << endl;
            cout << "y_  = " << y_ << endl;
            ASSERT(0);
        }

    }
}

template <typename CCS, typename FS>
void
mtv(int m, int n, int max_nnz,
    const GeCCSMatrix<CCS> &A, const GeMatrix<FS> &A_)
{
    typedef typename GeCCSMatrix<CCS>::ElementType  ElementType;

    DenseVector<Array<ElementType> >  x(m), y, y_;

    for (int i=1; i<=m; ++i) {
        x(i) = rand() % 10;
    }

    y  = transpose(A)  * x;
    y_ = transpose(A_) * x;

    if (! lapack::isIdentical(y, y_, "y", "y_")) {
        cerr << endl << "failed: y = A*x" << endl;
        ASSERT(0);
    }

    y  += transpose(A)  * x;
    y_ += transpose(A_) * x;

    if (! lapack::isIdentical(y, y_, "y", "y_")) {
        cerr << endl << "failed: y += A*x" << endl;
        ASSERT(0);
    }

    y  -= transpose(A)  * x;
    y_ -= transpose(A_) * x;

    if (! lapack::isIdentical(y, y_, "y", "y_")) {
        cerr << endl << "failed: y -= A*x" << endl;
        ASSERT(0);
    }

    ElementType  alpha, beta;
    for (int test=1; test<=20; ++test) {

        alpha = std::pow(2, 5 - std::max(1, rand() % 10));
        beta  = std::pow(2, 5 - std::max(1, rand() % 10));

        //
        // Reset y (and y_) to some random vector
        //
        for (int j=1; j<=n; ++j) {
            y(j) = rand() % 1000;
        }
        y_ = y;

        y  = beta*y  + alpha * transpose(A)  * x;
        y_ = beta*y_ + alpha * transpose(A_) * x;

        if (! lapack::isIdentical(y, y_, "y", "y_")) {
            cerr << endl << "failed: y = beta*y + alpha*A*x" << endl;
            cout << "alpha = " << alpha << endl;
            cout << "beta  = " << beta << endl;
            cout << "A_  = " << A_ << endl;
            cout << "x  = " << x << endl;
            cout << "y_  = " << y_ << endl;
            ASSERT(0);
        }

    }
}

int
main()
{
    srand(SEED);

    for (int run=1; run<=3000; ++run) {
        int m       = std::max(1, rand() % (MAX_M));
        int n       = std::max(1, rand() % (MAX_N));
        // check case 'nnz==0' at least onece
        int max_nnz = (run==1) ? 0 : (rand() % (MAX_NNZ));

        cerr << "run " << run << ":" << endl;

        for (int indexBase=-3; indexBase<=3; ++indexBase) {
            //
            //  Test non-square matrices
            //
            {
                cerr << "indexBase = " << indexBase << endl;
                cerr << "m =         " << m << endl;
                cerr << "n =         " << n << endl;
                cerr << "max_nnz =   " << max_nnz << endl << endl;

                GeCCSMatrix<CCS<double> >       A;
                GeMatrix<FullStorage<double> >  A_;

                setup(m, n, max_nnz, indexBase, A, A_);
                mv(m, n, max_nnz, A, A_);
                mtv(m, n, max_nnz, A, A_);
            }

            //
            //  Test square matrices
            //
            {
                cerr << "indexBase = " << indexBase << endl;
                cerr << "m x m = " << m << " x " << m << endl;
                cerr << "max_nnz = " << max_nnz << endl << endl;

                GeCCSMatrix<CCS<double> >       A;
                GeMatrix<FullStorage<double> >  A_;

                setup(m, m, max_nnz, indexBase, A, A_);
                mv(m, m, max_nnz, A, A_);
                mtv(m, m, max_nnz, A, A_);
            }
        }
    }
}
