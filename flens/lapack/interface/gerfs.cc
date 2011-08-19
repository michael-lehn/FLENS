#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define GERFS         sgerfs_
#   define GERFS_REF     sgerfs
#   define GERFS_NAME    "SGERFS"
#elif DOUBLE
#   define GERFS         dgerfs_
#   define GERFS_REF     dgerfs
#   define GERFS_NAME    "DGERFS"
#elif COMPLEX_SINGLE
#   define GERFS         cgerfs_
#   define GERFS_REF     cgerfs
#   define GERFS_NAME    "CGERFS"
#elif COMPLEX_DOUBLE
#   define GERFS         zgerfs_
#   define GERFS_REF     zgerfs
#   define GERFS_NAME    "ZGERFS"
#endif

using namespace flens;

extern "C" {

void
GERFS_REF(char *TRANS, INT *N, INT *NRHS, FLOAT *A, INT *LDA,
          FLOAT *AF, INT *LDAF, INT *IPIV, FLOAT *B, INT *LDB,
          FLOAT *X, INT *LDX, FLOAT *FERR, FLOAT *BERR,
          FLOAT *WORK, INT *IWORK, INT *INFO);


void
gerfsErrorCheck(char *TRANS, INT *N, INT *NRHS, FLOAT *A, INT *LDA,
                FLOAT *AF, INT *LDAF, INT *IPIV, FLOAT *B, INT *LDB,
                FLOAT *X, INT *LDX, FLOAT *FERR, FLOAT *BERR,
                FLOAT *WORK, INT *IWORK, INT *INFO)
{
    if ((*TRANS!='N') && (*TRANS!='n') && (*TRANS!='T') && (*TRANS!='t')
     && (*TRANS!='C') && (*TRANS!='c') && (*TRANS!='R') && (*TRANS!='r'))
    {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*NRHS<0) {
        *INFO = -3;
    } else if (*LDA<max(1,*N)) {
        *INFO = -5;
    } else if (*LDAF<max(1,*N)) {
        *INFO = -7;
    } else if (*LDB<max(1,*N)) {
        *INFO = -10;
    } else if (*LDX<max(1,*N)) {
        *INFO = -12;
    }
}


void
GERFS(char *TRANS, INT *N, INT *NRHS, FLOAT *A, INT *LDA,
      FLOAT *AF, INT *LDAF, INT *IPIV, FLOAT *B, INT *LDB,
      FLOAT *X, INT *LDX, FLOAT *FERR, FLOAT *BERR,
      FLOAT *WORK, INT *IWORK, INT *INFO)
{
    gerfsErrorCheck(TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,
                    X, LDX, FERR, BERR, WORK, IWORK, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(GERFS_NAME, INFO, strlen(GERFS_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifndef USE_GERFS_REF
    ASSERT(0);  // not yet implemented
#   else
    GERFS_REF(TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,
              X, LDX, FERR, BERR, WORK, IWORK, INFO);
#   endif
}

} // extern "C"

