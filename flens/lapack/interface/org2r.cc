#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define ORG2R         sorg2r_
#   define ORG2R_REF     sorg2r
#   define ORG2R_NAME    "SORG2R"
#elif DOUBLE
#   define ORG2R         dorg2r_
#   define ORG2R_REF     dorg2r
#   define ORG2R_NAME    "DORG2R"
#elif COMPLEX_SINGLE
#   define ORG2R         corg2r_
#   define ORG2R_REF     corg2r
#   define ORG2R_NAME    "CORG2R"
#elif COMPLEX_DOUBLE
#   define ORG2R         zorg2r_
#   define ORG2R_REF     zorg2r
#   define ORG2R_NAME    "ZORG2R"
#endif

using namespace flens;

extern "C" {

void
ORG2R_REF(INT *M, INT *N, INT *K, FLOAT *A, INT *LDA, FLOAT *TAU,
          FLOAT *WORK, INT *LWORK, INT *INFO);

void
org2rErrorCheck(INT *M, INT *N, INT *K, FLOAT *A, INT *LDA, FLOAT *TAU,
                FLOAT *WORK, INT *LWORK, INT *INFO)
{
    bool lQuery = (*LWORK==-1);
    
    if (*M<0) {
        *INFO = -1;
    } else if ((*N<0) || (*N>*M)) {
        *INFO = -2;
    } else if ((*K<0) || (*K>*N)) {
        *INFO = -3;
    } else if (*LDA<std::max(INT(1), *M)) {
        *INFO = -5;
    } else if ((*LWORK<std::max(INT(1), *N)) && (!lQuery)) {
        *INFO = -8;
    }
}


void
ORG2R(INT *M, INT *N, INT *K, FLOAT *A, INT *LDA, FLOAT *TAU,
      FLOAT *WORK, INT *LWORK, INT *INFO)
{
#   ifdef DEBUG_INTERFACE
    std::cerr << ORG2R_NAME << std::endl;
#   endif

    org2rErrorCheck(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(ORG2R_NAME, INFO, strlen(ORG2R_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifndef USE_ORG2R_REF
    ASSERT(0);  // not yet implemented
#   else
    ORG2R_REF(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);
#   endif
}

} // extern "C"

