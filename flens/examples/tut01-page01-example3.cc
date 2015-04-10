#include <cstdio>
#include <cstdlib>

struct DGeMatrix
{
    double *data;
    long   numRows, numCols;
};

struct DGeMatrixView
{
    double *data;
    long   numRows, numCols;
    long   leadingDimension;
};

// Constructor for DGeMatrix
void
DGeMatrix_alloc(DGeMatrix *A, long m, long n)
{
    A->data    = (double *)malloc(m*n*sizeof(double));
    A->numRows = m;
    A->numCols = n;
}

// Destructor for DGeMatrix
void
DGeMatrix_free(DGeMatrix *A)
{
    free(A->data);
    A->data = 0;
}

// View from DGematrix
void
DGeMatrix_view(const DGeMatrix *A, DGeMatrixView *B,
               long fromRow, long toRow,
               long fromCol, long toCol)
{
    long ldA = A->numRows;

    B->data             = A->data + fromRow + fromCol*ldA;
    B->numRows          = toRow-fromRow+1;
    B->numCols          = toCol-fromCol+1;
    B->leadingDimension = ldA;
}

// Change entry in DGeMatrix
void
DGeMatrix_set(const DGeMatrix *A, long row, long col, double value)
{
    A->data[row+col*A->numRows] = value;
}

// Get entry from DGeMatrix
double
DGeMatrix_get(const DGeMatrix *A, long row, long col)
{
    return A->data[row+col*A->numRows];
}

// Change entry in DGeMatrixView
void
DGeMatrixView_set(const DGeMatrixView *A, long row, long col, double value)
{
    A->data[row+col*A->leadingDimension] = value;
}

// Get entry from DGeMatrixView
double
DGeMatrixView_get(const DGeMatrixView *A, long row, long col)
{
    return A->data[row+col*A->leadingDimension];
}

int
main()
{
    // DGeMatrix  A(3,4);
    DGeMatrix A;
    DGeMatrix_alloc(&A, 3, 4);

    // A = 1,  2,  3,  4,
    //     5,  6,  7,  8,
    //     9, 10, 11, 12;
    //
    // The list initializer always fills the matrix row wise.  Even if col wise
    // would be faster.  But anyway, the sole purpose of the list initializer is
    // for writing compact examples in the tutorial
    //
    long counter = 1;
    for (long j=0; j<A.numCols; ++j) {
        for (long i=0; i<A.numRows; ++i) {
            DGeMatrix_set(&A, i, j, counter);
            ++counter;
        }
    }

    // cout << "A = " << A << endl;
    printf("A =\n");
    for (long i=0; i<A.numRows; ++i) {
        for (long j=0; j<A.numCols; ++j) {
            printf("%8.3lf", DGeMatrix_get(&A, i, j));
        }
        printf("\n");
    }
    printf("\n");


    // DGeMatrixView B = A(_(2,3),_);
    DGeMatrixView B;
    DGeMatrix_view(&A, &B, 1, 2, 0, 3);

    // cout << "B = " << B << endl;
    printf("B =\n");
    for (long i=0; i<B.numRows; ++i) {
        for (long j=0; j<B.numCols; ++j) {
            printf("%8.3lf", DGeMatrixView_get(&B, i, j));
        }
        printf("\n");
    }
    printf("\n");

    // B(1,2) = 42;
    DGeMatrixView_set(&B, 0, 1, 42);

    // cout << "Changed B(1,2)." << std::endl;
    printf("Changed B(1,2).\n");

    // cout << "A = " << A << endl;
    printf("A =\n");
    for (long i=0; i<A.numRows; ++i) {
        for (long j=0; j<A.numCols; ++j) {
            printf("%8.3lf", DGeMatrix_get(&A, i, j));
        }
        printf("\n");
    }
    printf("\n");

    // cout << "B = " << B << endl;
    printf("B =\n");
    for (long i=0; i<B.numRows; ++i) {
        for (long j=0; j<B.numCols; ++j) {
            printf("%8.3lf", DGeMatrixView_get(&B, i, j));
        }
        printf("\n");
    }
    printf("\n");

    // Destructor for A
    DGeMatrix_free(&A);
}
