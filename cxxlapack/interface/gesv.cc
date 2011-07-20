#include <cxxlapack/cxxlapack.cxx>
#include <cxxblas/cxxblas.cxx>

#include <cxxlapack/interface/aux.h>

#ifdef SINGLE
#   define GECON        sgecon_
#   define GEEQU        sgeequ_
#   define GERFS        sgerfs_
#   define GESV         sgesv_
#   define GESVX        sgesvx_
#   define GETF2        sgetf2_
#   define GETRF        sgetrf_
#   define GETRI        sgetri_
#   define GETRS        sgetrs_
#elif DOUBLE
#   define GETRF        dgetrf_
#elif COMPLEX_SINGLE
#   define GETRF        cgetrf_
#elif COMPLEX_DOUBLE
#   define GETRF        zgetrf_
#endif

#ifdef SINGLE
#   define GECON_REF    sgecon
#   define GEEQU_REF    sgeequ
#   define GERFS_REF    sgerfs
#   define GESV_REF     sgesv
#   define GESVX_REF    sgesvx
#   define GETF2_REF    sgetf2
#   define GETRF_REF    sgetrf
#   define GETRI_REF    sgetri
#   define GETRS_REF    sgetrs
#endif

using cxxlapack::gecon;
using cxxlapack::geequ;
using cxxlapack::gerfs;
using cxxlapack::gesv;
//using cxxlapack::gesvx;
using cxxlapack::getf2;
using cxxlapack::getrf;
using cxxlapack::getri;
using cxxlapack::getrs;

extern "C" {

//- gecon ----------------------------------------------------------------------

#ifdef USE_GECON_REF
void
GECON_REF(const char *NORM, const INT *N, const FLOAT *A, const INT *LDA,
          const FLOAT *ANORM, FLOAT *RCOND, FLOAT *WORK, INT *IWORK,
          INT *INFO);
#endif

void
GECON(const char *NORM, const INT *N, const FLOAT *A, const INT *LDA,
      const FLOAT *ANORM, FLOAT *RCOND, FLOAT *WORK, INT *IWORK,
      INT *INFO)
{
    bool checkNorm = false;
    cxxlapack::Norm norm;

    if ((*NORM=='1') || (*NORM=='O')) {
        norm = cxxlapack::OneNorm;
        checkNorm = true;
    }
    if (*NORM=='I') {
        norm = cxxlapack::InfinityNorm;
        checkNorm = true;
    }

    if (checkNorm) {
        *INFO = gecon(ColMajor, norm, *N, A, *LDA,
                      *ANORM, *RCOND, WORK, IWORK);
    } else {
        *INFO = -1;
    }
    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA("SGECON", INFO, 6);
        return;
    }

#   ifdef USE_GECON_REF
    GECON_REF(NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, INFO);
#   endif
}

//- geequ ----------------------------------------------------------------------

#ifdef USE_GEEQU_REF
void
GEEQU_REF(const INT *M, const INT *N, const FLOAT *A, const INT *LDA,
          FLOAT *R, FLOAT *C, FLOAT *ROWCND, FLOAT *COLCND, FLOAT *AMAX,
          INT *INFO);
#endif

void
GEEQU(const INT *M, const INT *N, const FLOAT *A, const INT *LDA,
      FLOAT *R, FLOAT *C, FLOAT *ROWCND, FLOAT *COLCND, FLOAT *AMAX,
      INT *INFO)
{
    *INFO = geequ(ColMajor, *M, *N, A, *LDA, R, C, *ROWCND, *COLCND, *AMAX);
    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA("SGEEQU", INFO, 6);
        return;
    }
#   ifdef USE_GEEQU_REF
    GEEQU_REF(M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO);
#   endif
}

//- gerfs ----------------------------------------------------------------------

#ifdef USE_GERFS_REF
void
GERFS_REF(const char *TRANS, const INT *N, const INT *NRHS,
          const FLOAT *A, const INT *LDA,
          const FLOAT *AF, const INT *LDAF,
          const INT *IPIV,
          const FLOAT *B, const INT *LDB,
          FLOAT *X, INT *LDX,
          FLOAT *FERR, FLOAT *BERR,
          FLOAT *WORK, INT *IWORK,
          INT *INFO);
#endif

void
GERFS(const char *TRANS, const INT *N, const INT *NRHS,
      const FLOAT *A, const INT *LDA,
      const FLOAT *AF, const INT *LDAF,
      const INT *IPIV,
      const FLOAT *B, const INT *LDB,
      FLOAT *X, INT *LDX,
      FLOAT *FERR, FLOAT *BERR,
      FLOAT *WORK, INT *IWORK,
      INT *INFO)
{
    bool checkTransA = false;
    Transpose transA = NoTrans;

    if ((*TRANS=='N') || (*TRANS=='n')) {
        checkTransA = true;
    }
    if ((*TRANS=='T') || (*TRANS=='t')) {
        transA = Trans;
        checkTransA = true;
    }
    if ((*TRANS=='C') || (*TRANS=='c')) {
        transA = ConjTrans;
        checkTransA = true;
    }
    if ((*TRANS=='R') || (*TRANS=='r')) {
        transA = Conj;
        checkTransA = true;
    }

    if (checkTransA) {
        *INFO = gerfs(ColMajor,
                      transA, *N, *NRHS,
                      reinterpret_cast<const CXXFLOAT *>(A), *LDA,
                      reinterpret_cast<const CXXFLOAT *>(AF), *LDAF,
                      IPIV,
                      reinterpret_cast<const CXXFLOAT *>(B), *LDB,
                      reinterpret_cast<CXXFLOAT *>(X), *LDX,
                      reinterpret_cast<CXXFLOAT *>(FERR),
                      reinterpret_cast<CXXFLOAT *>(BERR),
                      reinterpret_cast<CXXFLOAT *>(WORK), IWORK);
    } else {
        *INFO = -1;
    }

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA("SGERFS", INFO, 6);
        return;
    }
#   ifdef USE_GERFS_REF
    GERFS_REF(TRANS, N, NRHS, A, LDA, AF, LDAF,
              IPIV, B, LDB, X, LDX, FERR, BERR,
              WORK, IWORK, INFO);
#   endif
}

//- gesv -----------------------------------------------------------------------

#ifdef USE_GESV_REF
void
GESV_REF(const INT *N, const INT *NRHS,
         FLOAT *A, const INT *LDA,
         INT *IPIV,
         FLOAT *B, const INT *LDB,
         INT *INFO);
#endif

void
GESV(const INT *N, const INT *NRHS,
     FLOAT *A, const INT *LDA,
     INT *IPIV,
     FLOAT *B, const INT *LDB,
     INT *INFO)
{
    *INFO = gesv(ColMajor, *N, *NRHS, A, *LDA, IPIV, B, *LDB);
    for (INT i=0; i<*N; ++i) {
        ++IPIV[i];
    }
    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA("SGESV", INFO, 5);
        *INFO = -(*INFO);
    }
    if (*INFO<0) {
        return;
    }
#   ifdef USE_GESV_REF
    GESV_REF(N, NRHS, A, LDA, IPIV, B, LDB, INFO);
#   endif
}

//- gesvx ----------------------------------------------------------------------

#ifdef USE_GESVX_REF
void
GESVX_REF(const char *FACT, const char *TRANS, const INT *N, const INT *NRHS,
          FLOAT *A, const INT *LDA, FLOAT *AF, const INT *LDAF,
          INT *IPIV,
          char *EQUED,
          FLOAT *R, FLOAT *C,
          FLOAT *B, const INT *LDB,
          FLOAT *X, const INT *LDX,
          FLOAT *RCOND,
          FLOAT *FERR, FLOAT *BERR,
          FLOAT *WORK, INT *IWORK, INT *INFO);
#endif

void
GESVX(const char *FACT, const char *TRANS, const INT *N, const INT *NRHS,
      FLOAT *A, const INT *LDA, FLOAT *AF, const INT *LDAF,
      INT *IPIV,
      char *EQUED,
      FLOAT *R, FLOAT *C,
      FLOAT *B, const INT *LDB,
      FLOAT *X, const INT *LDX,
      FLOAT *RCOND,
      FLOAT *FERR, FLOAT *BERR,
      FLOAT *WORK, INT *IWORK, INT *INFO)
{
//
//  Test the input parameters.
//
    FLOAT ONE = 1;
    FLOAT SMLNUM, BIGNUM;

    *INFO = 0;
    bool NOFACT = ((*FACT=='N') || (*FACT=='n'));
    bool EQUIL  = ((*FACT=='E') || (*FACT=='e'));
    bool NOTRAN = ((*TRANS=='N') || (*TRANS=='n'));
    bool ROWEQU;
    bool COLEQU;
    if (NOFACT || EQUIL) {
        *EQUED = 'N';
        ROWEQU = false;
        COLEQU = false;
    } else {
        ROWEQU = (*EQUED=='R') || (*EQUED=='r')
              || (*EQUED=='B') || (*EQUED=='b');
        COLEQU = (*EQUED=='C') || (*EQUED=='c')
              || (*EQUED=='B') || (*EQUED=='b');
        SMLNUM = cxxlapack::lamch<FLOAT>(cxxlapack::SafeMin);
        BIGNUM = ONE / SMLNUM;
    }

    if (!NOFACT && !EQUIL && (toupper(*FACT)!='F') ) {
        *INFO = -1;
    } else if (!NOTRAN && (toupper(*TRANS)!='T') && (toupper(*TRANS)!='C')) {
        *INFO = -2;
    } else if (*N<0) {
        *INFO = -3;
    } else if (*NRHS<0) {
        *INFO = -4;
    } else if (*LDA<max(1,*N)) {
        *INFO = -6;
    } else if (*LDAF<max(1,*N)) {
        *INFO = -8;
    } else if ((toupper(*FACT)=='F')
           && !(ROWEQU || COLEQU || toupper(*EQUED)=='N'))
    {
        *INFO = -10;
    } else {
         if (ROWEQU) {
            FLOAT RCMIN = BIGNUM;
            FLOAT RCMAX = FLOAT(0);
            for (INT j=0; j<*N; ++j) {
                RCMIN = min(RCMIN, R[j]);
                RCMAX = max(RCMAX, R[j]);
            }
            if (RCMIN<=FLOAT(0)) {
                *INFO = -11;
            } else if (*N>0) {
            //    ROWCND = max(RCMIN, SMLNUM) / min(RCMAX, BIGNUM);
            } else {
            //    ROWCND = FLOAT(1);
            }
        }
        if (COLEQU && (*INFO==0)) {
            FLOAT RCMIN = BIGNUM;
            FLOAT RCMAX = FLOAT(0);
            for (INT j=0; j<*N; ++j) {
                RCMIN = min(RCMIN, C[j]);
                RCMAX = max(RCMAX, C[j]);
            }
            if (RCMIN<=FLOAT(0)) {
                *INFO = -12;
            } else if (*N>0) {
            //    COLCND = max(RCMIN, SMLNUM) / min(RCMAX, BIGNUM);
            } else {
            //    COLCND = FLOAT(1);
            }
            if (*INFO==0) {
                if (*LDB<max(1, *N)) {
                    *INFO = -14;
                } else if (*LDX<max(1,*N)) {
                    *INFO = -16;
                }
            }
        }
    }

    if (*LDB<max(1, *N)) {
        *INFO = -14;
    }
    if (*LDX<max(1,*N)) {
        *INFO = -16;
    }

    if (*INFO!=0) {
        *INFO = -(*INFO);
        XERBLA("SGESVX", INFO, 6);
        return;
    }
    GESVX_REF(FACT, TRANS, N, NRHS, A, LDA, AF, LDAF,
              IPIV, EQUED, R, C, B, LDB, X, LDX,
              RCOND, FERR, BERR, WORK, IWORK, INFO);
}

//- getf2 ----------------------------------------------------------------------

#ifdef USE_GETF2_REF
void
GETF2_REF(const INT *M, const INT *N,
          FLOAT *A, const INT *LDA, INT *IPIV,
          INT *INFO);
#endif

void
GETF2(const INT *M, const INT *N,
      FLOAT *A, const INT *LDA, INT *IPIV,
      INT *INFO)
{
    *INFO = getf2(ColMajor,
                  *M, *N,
                  reinterpret_cast<CXXFLOAT *>(A), *LDA,
                  IPIV);
    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA("SGETF2", INFO, 6);
        *INFO = -(*INFO);
    }
    for (INT i=0; i<std::min(*M, *N); ++i) {
        ++IPIV[i];
    }
    if (*INFO<0) {
        return;
    }

#   ifdef USE_GETF2_REF
    GETF2_REF(M, N, A, LDA, IPIV, INFO);
#   endif

}

//- getrf ----------------------------------------------------------------------

#ifdef USE_GETRF_REF
void
GETRF_REF(const INT *M, const INT *N, FLOAT *A, const INT *LDA,
          INT *IPIV, INT *INFO);
#endif

void
GETRF(const INT *M, const INT *N, FLOAT *A, const INT *LDA,
      INT *IPIV, INT *INFO)
{
    *INFO = getrf(ColMajor,
                  *M, *N,
                  reinterpret_cast<CXXFLOAT *>(A), *LDA,
                  IPIV);
    
    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA("SGETRF", INFO, 6);
        *INFO = -(*INFO);
    }
    for (INT i=0; i<*M; ++i) {
        ++IPIV[i];
    }
    if (*INFO<0) {
        return;
    }
#   ifdef USE_GETRF_REF
    GETRF_REF(M, N, A, LDA, IPIV, INFO);
#   endif
}

//- getri ----------------------------------------------------------------------

#ifdef USE_GETRI_REF
void
GETRI_REF(const INT *N, FLOAT *A, const INT *LDA, const INT *IPIV,
          FLOAT *WORK, const INT *LWORK, INT *INFO);
#endif

void
GETRI(const INT *N, FLOAT *A, const INT *LDA, const INT *IPIV,
      FLOAT *WORK, const INT *LWORK, INT *INFO)
{
    INT *iPiv = const_cast<INT *>(IPIV);
    for (INT i=0; i<*N; ++i) {
        --iPiv[i];
    }
    *INFO = getri(ColMajor,
                  *N, reinterpret_cast<CXXFLOAT *>(A), *LDA,
                  IPIV,
                  reinterpret_cast<CXXFLOAT *>(WORK), *LWORK);
    
    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA("SGETRI", INFO, 6);
        *INFO = -(*INFO);
    }
    for (INT i=0; i<*N; ++i) {
        ++iPiv[i];
    }
    if (*INFO<0) {
        return;
    }
#   ifdef USE_GETRI_REF
    GETRI_REF(N, A, LDA, IPIV, WORK, LWORK, INFO);
#   endif
}

//- getrs ----------------------------------------------------------------------

#ifdef USE_GETRS_REF
void
GETRS_REF(const char *TRANS, const INT *N, const INT *NRHS,
          const FLOAT *A, const INT *LDA, const INT *IPIV,
          FLOAT *B, const INT *LDB, INT *INFO);
#endif

void
GETRS(const char *TRANS, const INT *N, const INT *NRHS,
      const FLOAT *A, const INT *LDA, const INT *IPIV,
      FLOAT *B, const INT *LDB, INT *INFO)
{
    bool checkTransA = false;
    Transpose transA = NoTrans;

    if ((*TRANS=='N') || (*TRANS=='n')) {
        checkTransA = true;
    }
    if ((*TRANS=='T') || (*TRANS=='t')) {
        transA = Trans;
        checkTransA = true;
    }
    if ((*TRANS=='C') || (*TRANS=='c')) {
        transA = ConjTrans;
        checkTransA = true;
    }
    if ((*TRANS=='R') || (*TRANS=='r')) {
        transA = Conj;
        checkTransA = true;
    }
//
//  change pivot indices to zero based indices
//
    INT *iPiv = const_cast<INT *>(IPIV);
    for (INT i=0; i<*N; ++i) {
        --iPiv[i];
    }
    if (checkTransA) {
        *INFO = getrs(ColMajor,
                      transA, *N, *NRHS,
                      reinterpret_cast<const CXXFLOAT *>(A), *LDA, IPIV,
                      reinterpret_cast<CXXFLOAT *>(B), *LDB);
    } else {
        *INFO = -1;
    }
    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA("SGETRS", INFO, 6);
        *INFO = -(*INFO);
    }
//
//  restore pivot indices
//
    for (INT i=0; i<*N; ++i) {
        ++iPiv[i];
    }
    if (*INFO<0) {
        return;
    }
#ifdef USE_GETRS_REF
    GETRS_REF(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO);
#endif
}

// dummy

#ifndef DUMMY

// void slamch_() { fprintf(stderr, "not implemented1\n"); }
// void sgetrs_() { fprintf(stderr, "not implemented2\n"); }
void spotrf_() { fprintf(stderr, "not implemented3\n"); }
void spotri_() { fprintf(stderr, "not implemented4\n"); }
void strtri_() { fprintf(stderr, "not implemented5\n"); }
void strti2_() { fprintf(stderr, "not implemented6\n"); }
// void slaswp_() { fprintf(stderr, "not implemented7\n"); }
// void sgesv_()  { fprintf(stderr, "not implemented8\n"); }
void spotf2_() { fprintf(stderr, "not implemented9\n"); }




void sgbtrf_() { fprintf(stderr, "not implemented10\n"); }


void ieeeck_() { fprintf(stderr, "not implemented11\n"); }
// void second_() { fprintf(stderr, "not implemented12\n"); }
// void ilaver_() { fprintf(stderr, "not implemented13\n"); }
// void sgeequ_() { fprintf(stderr, "not implemented14\n"); }

void sgbequ_() {
#   ifdef SHOW_TODO
    fprintf(stderr, "not implemented15\n");
#   endif
}

void spoequ_() {
#   ifdef SHOW_TODO
    fprintf(stderr, "not implemented16\n");
#   endif
}

void sppequ_() {
#   ifdef SHOW_TODO
    fprintf(stderr, "not implemented17\n");
#   endif
}

void spbequ_() {
#   ifdef SHOW_TODO
    fprintf(stderr, "not implemented18\n");
#   endif
}

//  void slacpy_() { fprintf(stderr, "not implemented19\n"); }
void slangb_() { fprintf(stderr, "not implemented20\n"); }
// void slaset_() { fprintf(stderr, "not implemented21\n"); }
void sgbtrs_() { fprintf(stderr, "not implemented22\n"); }
// void slange_() { fprintf(stderr, "not implemented23\n"); }
void sgbrfs_() { fprintf(stderr, "not implemented24\n"); }
void sgbcon_() { fprintf(stderr, "not implemented25\n"); }
// void sgetri_() { fprintf(stderr, "not implemented26\n"); }
// void sgerfs_() { fprintf(stderr, "not implemented27\n"); }
// void sgecon_() { fprintf(stderr, "not implemented28\n"); }
// void slarnv_() { fprintf(stderr, "not implemented29\n"); }
void sgttrf_() { fprintf(stderr, "not implemented30\n"); }
void slangt_() { fprintf(stderr, "not implemented31\n"); }
void sgttrs_() { fprintf(stderr, "not implemented32\n"); }
void sgtcon_() { fprintf(stderr, "not implemented33\n"); }
void slagtm_() { fprintf(stderr, "not implemented34\n"); }
void sgtrfs_() { fprintf(stderr, "not implemented35\n"); }
void spbtrf_() { fprintf(stderr, "not implemented36\n"); }
void spbtrs_() { fprintf(stderr, "not implemented37\n"); }
void slansb_() { fprintf(stderr, "not implemented38\n"); }
void spbrfs_() { fprintf(stderr, "not implemented39\n"); }
void spbcon_() { fprintf(stderr, "not implemented40\n"); }
void spotrs_() { fprintf(stderr, "not implemented41\n"); }
void sporfs_() { fprintf(stderr, "not implemented42\n"); }
void slansy_() { fprintf(stderr, "not implemented43\n"); }
void spocon_() { fprintf(stderr, "not implemented44\n"); }
void spstrf_() { fprintf(stderr, "not implemented45\n"); }
void spptrf_() { fprintf(stderr, "not implemented46\n"); }
void spptri_() { fprintf(stderr, "not implemented47\n"); }
void spptrs_() { fprintf(stderr, "not implemented48\n"); }
void spprfs_() { fprintf(stderr, "not implemented49\n"); }
void slansp_() { fprintf(stderr, "not implemented50\n"); }
void sppcon_() { fprintf(stderr, "not implemented51\n"); }
void spttrf_() { fprintf(stderr, "not implemented52\n"); }
void slanst_() { fprintf(stderr, "not implemented53\n"); }
void spttrs_() { fprintf(stderr, "not implemented54\n"); }
void sptrfs_() { fprintf(stderr, "not implemented55\n"); }
void sptcon_() { fprintf(stderr, "not implemented56\n"); }
void sgeqp3_() { fprintf(stderr, "not implemented57\n"); }
void sgeqpf_() { fprintf(stderr, "not implemented58\n"); }
void ssptrf_() { fprintf(stderr, "not implemented59\n"); }
void ssptri_() { fprintf(stderr, "not implemented60\n"); }
void ssptrs_() { fprintf(stderr, "not implemented61\n"); }
void ssprfs_() { fprintf(stderr, "not implemented62\n"); }
void sspcon_() { fprintf(stderr, "not implemented63\n"); }
void ssytrf_() { fprintf(stderr, "not implemented64\n"); }
void ssytri2_() { fprintf(stderr, "not implemente65d\n"); }
void ssytrs_() { fprintf(stderr, "not implemented66\n"); }
void ssytrs2_() { fprintf(stderr, "not implemented67\n"); }
void ssyrfs_() { fprintf(stderr, "not implemented68\n"); }
void ssycon_() { fprintf(stderr, "not implemented69\n"); }
void slantb_() { fprintf(stderr, "not implemented70\n"); }
// void slantr_() { fprintf(stderr, "not implemented71\n"); }
void stbtrs_() { fprintf(stderr, "not implemented72\n"); }
void stbrfs_() { fprintf(stderr, "not implemented73\n"); }
void stbcon_() { fprintf(stderr, "not implemented74\n"); }
void slatbs_() { fprintf(stderr, "not implemented75\n"); }
void stptri_() { fprintf(stderr, "not implemented76\n"); }
void slantp_() { fprintf(stderr, "not implemented77\n"); }
void stptrs_() { fprintf(stderr, "not implemented78\n"); }
void stprfs_() { fprintf(stderr, "not implemented79\n"); }
void stpcon_() { fprintf(stderr, "not implemented80\n"); }
void slatps_() { fprintf(stderr, "not implemented81\n"); }
void strtrs_() { fprintf(stderr, "not implemented82\n"); }
void strrfs_() { fprintf(stderr, "not implemented83\n"); }
void strcon_() { fprintf(stderr, "not implemented84\n"); }
void slatrs_() { fprintf(stderr, "not implemented85\n"); }
void sgeqr2_() { fprintf(stderr, "not implemented86\n"); }
void stzrqf_() { fprintf(stderr, "not implemented87\n"); }
void stzrzf_() { fprintf(stderr, "not implemented88\n"); }
void sgtsv_() { fprintf(stderr, "not implemented89\n"); }
void sgtsvx_() { fprintf(stderr, "not implemented90\n"); }
void sgels_() { fprintf(stderr, "not implemented91\n"); }
void sgelsx_() { fprintf(stderr, "not implemented92\n"); }
void sgelsy_() { fprintf(stderr, "not implemented93\n"); }
void sgelss_() { fprintf(stderr, "not implemented94\n"); }
void sgelsd_() { fprintf(stderr, "not implemented95\n"); }
void slaqsb_() { fprintf(stderr, "not implemented96\n"); }
void spbsv_() { fprintf(stderr, "not implemented97\n"); }
void spbsvx_() { fprintf(stderr, "not implemented98\n"); }
void slaqsp_() { fprintf(stderr, "not implemented99\n"); }
void sppsv_() { fprintf(stderr, "not implemented100\n"); }
void sppsvx_() { fprintf(stderr, "not implemented101\n"); }
void sptsv_() { fprintf(stderr, "not implemented102\n"); }
void sptsvx_() { fprintf(stderr, "not implemented103\n"); }
void sspsv_() { fprintf(stderr, "not implemented104\n"); }
void sspsvx_() { fprintf(stderr, "not implemented105\n"); }
void sgelqf_() { fprintf(stderr, "not implemented106\n"); }
void sgelq2_() { fprintf(stderr, "not implemented107\n"); }
void sorglq_() { fprintf(stderr, "not implemented108\n"); }
void sorgl2_() { fprintf(stderr, "not implemented109\n"); }
void sormlq_() { fprintf(stderr, "not implemented110\n"); }
void sorml2_() { fprintf(stderr, "not implemented111\n"); }
void spstf2_() { fprintf(stderr, "not implemented112\n"); }
void sgeqlf_() { fprintf(stderr, "not implemented113\n"); }
void sgeql2_() { fprintf(stderr, "not implemented114\n"); }
void sorgql_() { fprintf(stderr, "not implemented115\n"); }
void sorg2l_() { fprintf(stderr, "not implemented116\n"); }
void sormql_() { fprintf(stderr, "not implemented117\n"); }
void sorm2l_() { fprintf(stderr, "not implemented118\n"); }
void sgeqrf_() { fprintf(stderr, "not implemented119\n"); }
void sgeqrfp_() { fprintf(stderr, "not implemented120\n"); }
void sgeqr2p_() { fprintf(stderr, "not implemented121\n"); }
void sorgqr_() { fprintf(stderr, "not implemented122\n"); }
void sorg2r_() { fprintf(stderr, "not implemented123\n"); }
void sormqr_() { fprintf(stderr, "not implemented124\n"); }
void sorm2r_() { fprintf(stderr, "not implemented125\n"); }
void sgerqf_() { fprintf(stderr, "not implemented126\n"); }
void sgerq2_() { fprintf(stderr, "not implemented127\n"); }
void sorgrq_() { fprintf(stderr, "not implemented128\n"); }
void sorgr2_() { fprintf(stderr, "not implemented129\n"); }
void sormrq_() { fprintf(stderr, "not implemented130\n"); }
void sormr2_() { fprintf(stderr, "not implemented131\n"); }
void slanhs_() { fprintf(stderr, "not implemented132\n"); }
// void slabad_() { fprintf(stderr, "not implemented133\n"); }
void slascl_() { fprintf(stderr, "not implemented134\n"); }
void sgebd2_() { fprintf(stderr, "not implemented135\n"); }
void sbdsqr_() { fprintf(stderr, "not implemented136\n"); }
void slarf_() { fprintf(stderr, "not implemented137\n"); }
void sormrz_() { fprintf(stderr, "not implemented138\n"); }
void slatzm_() { fprintf(stderr, "not implemented139\n"); }
//void sgesvx_() { fprintf(stderr, "not implemented140\n"); }
void sgbsv_() { fprintf(stderr, "not implemented141\n"); }
void sgbsvx_() { fprintf(stderr, "not implemented142\n"); }
void sposv_() { fprintf(stderr, "not implemented143\n"); }
void sposvx_() { fprintf(stderr, "not implemented144\n"); }
void ssysv_() { fprintf(stderr, "not implemented145\n"); }
void ssysvx_() { fprintf(stderr, "not implemented146\n"); }
// void slaqge_() { fprintf(stderr, "not implemented147\n"); }
void sgbtf2_() { fprintf(stderr, "not implemented148\n"); }
void slaqgb_() { fprintf(stderr, "not implemented149\n"); }
void slaqsy_() { fprintf(stderr, "not implemented150\n"); }
void ssytf2_() { fprintf(stderr, "not implemented151\n"); }
void ssytri_() { fprintf(stderr, "not implemented152\n"); }
void spbtf2_() { fprintf(stderr, "not implemented153\n"); }
// void slartg_() { fprintf(stderr, "not implemented154\n"); }

void sgbbrd_() { fprintf(stderr, "not implemented155\n"); }
void sgebrd_() { fprintf(stderr, "not implemented156\n"); }
void sorgbr_() { fprintf(stderr, "not implemented157\n"); }
void sbdsdc_() { fprintf(stderr, "not implemented158\n"); }
void sgebak_() { fprintf(stderr, "not implemented159\n"); }
void sgebal_() { fprintf(stderr, "not implemented160\n"); }
void slarfg_() { fprintf(stderr, "not implemented161\n"); }
void sgghrd_() { fprintf(stderr, "not implemented162\n"); }
void shgeqz_() { fprintf(stderr, "not implemented163\n"); }
void stgevc_() { fprintf(stderr, "not implemented164\n"); }
void sggbal_() { fprintf(stderr, "not implemented165\n"); }
void sggbak_() { fprintf(stderr, "not implemented166\n"); }
void sgehrd_() { fprintf(stderr, "not implemented167\n"); }
void sorghr_() { fprintf(stderr, "not implemented168\n"); }
void shseqr_() { fprintf(stderr, "not implemented169\n"); }
void strevc_() { fprintf(stderr, "not implemented170\n"); }
void shsein_() { fprintf(stderr, "not implemented171\n"); }
void sormhr_() { fprintf(stderr, "not implemented172\n"); }
void ssbtrd_() { fprintf(stderr, "not implemented173\n"); }
void ssytrd_() { fprintf(stderr, "not implemented174\n"); }
void sorgtr_() { fprintf(stderr, "not implemented175\n"); }
void ssptrd_() { fprintf(stderr, "not implemented176\n"); }
void sopgtr_() { fprintf(stderr, "not implemented177\n"); }
void ssteqr_() { fprintf(stderr, "not implemented178\n"); }
void ssterf_() { fprintf(stderr, "not implemented179\n"); }
void spteqr_() { fprintf(stderr, "not implemented180\n"); }
void sstebz_() { fprintf(stderr, "not implemented181\n"); }
void sstein_() { fprintf(stderr, "not implemented182\n"); }
void sstedc_() { fprintf(stderr, "not implemented183\n"); }
void sstemr_() { fprintf(stderr, "not implemented184\n"); }
void sorcsd_() { fprintf(stderr, "not implemented185\n"); }
void sgges_() { fprintf(stderr, "not implemented186\n"); }
void sggev_() { fprintf(stderr, "not implemented187\n"); }
void sggesx_() { fprintf(stderr, "not implemented188\n"); }
void sgesvd_() { fprintf(stderr, "not implemented189\n"); }
void sggevx_() { fprintf(stderr, "not implemented190\n"); }
void sgesdd_() { fprintf(stderr, "not implemented191\n"); }
void sgesvj_() { fprintf(stderr, "not implemented192\n"); }
void sgejsv_() { fprintf(stderr, "not implemented193\n"); }
void sgees_() { fprintf(stderr, "not implemented194\n"); }
void sgeev_() { fprintf(stderr, "not implemented195\n"); }
void slapy2_() { fprintf(stderr, "not implemented196\n"); }
void sgegs_() { fprintf(stderr, "not implemented197\n"); }
void sgegv_() { fprintf(stderr, "not implemented198\n"); }
void ssygv_() { fprintf(stderr, "not implemented199\n"); }
void ssygvd_() { fprintf(stderr, "not implemented200\n"); }
void ssygvx_() { fprintf(stderr, "not implemented201\n"); }
void sspgv_() { fprintf(stderr, "not implemented202\n"); }
void sspgvd_() { fprintf(stderr, "not implemented203\n"); }
void sspgvx_() { fprintf(stderr, "not implemented204\n"); }
void ssbgv_() { fprintf(stderr, "not implemented205\n"); }
void ssbgvd_() { fprintf(stderr, "not implemented206\n"); }
void ssbgvx_() { fprintf(stderr, "not implemented207\n"); }
void sstev_() { fprintf(stderr, "not implemented208\n"); }
void sstevx_() { fprintf(stderr, "not implemented209\n"); }
void sstevr_() { fprintf(stderr, "not implemented210\n"); }
void sstevd_() { fprintf(stderr, "not implemented211\n"); }
void ssyev_() { fprintf(stderr, "not implemented212\n"); }
void ssyevx_() { fprintf(stderr, "not implemented213\n"); }
void sspev_() { fprintf(stderr, "not implemented214\n"); }
void sspevx_() { fprintf(stderr, "not implemented215\n"); }
void ssbev_() { fprintf(stderr, "not implemented 216\n"); }
void ssbevx_() { fprintf(stderr, "not implemented 217\n"); }
void ssyevd_() { fprintf(stderr, "not implemented 218\n"); }
void sspevd_() { fprintf(stderr, "not implemented 219\n"); }
void ssbevd_() { fprintf(stderr, "not implemented 220\n"); }
void ssyevr_() { fprintf(stderr, "not implemented 221\n"); }
void sormbr_() { fprintf(stderr, "not implemented 222\n"); }
void strsyl_() { fprintf(stderr, "not implemented 223\n"); }
void strexc_() { fprintf(stderr, "not implemented 224\n"); }
void strsna_() { fprintf(stderr, "not implemented 225\n"); }
void strsen_() { fprintf(stderr, "not implemented 226\n"); }
void sgeevx_() { fprintf(stderr, "not implemented 227\n"); }
void sgeesx_() { fprintf(stderr, "not implemented 228\n"); }
void sggsvd_() { fprintf(stderr, "not implemented 229\n"); }
void sggsvp_() { fprintf(stderr, "not implemented 230\n"); }
void stgsja_() { fprintf(stderr, "not implemented 231\n"); }
void sggglm_() { fprintf(stderr, "not implemented 232\n"); }
void sgglse_() { fprintf(stderr, "not implemented 233\n"); }
void sggqrf_() { fprintf(stderr, "not implemented 234\n"); }
void sggrqf_() { fprintf(stderr, "not implemented 235\n"); }
void stgexc_() { fprintf(stderr, "not implemented 236\n"); }
void stgsen_() { fprintf(stderr, "not implemented 237\n"); }
void stgsna_() { fprintf(stderr, "not implemented 238\n"); }
void stgsyl_() { fprintf(stderr, "not implemented 239\n"); }
void sormtr_() { fprintf(stderr, "not implemented 240\n"); }
void sopmtr_() { fprintf(stderr, "not implemented 241\n"); }
void slaln2_() { fprintf(stderr, "not implemented 242\n"); }
void slasy2_() { fprintf(stderr, "not implemented 243\n"); }
void slanv2_() { fprintf(stderr, "not implemented 244\n"); }
void slaexc_() { fprintf(stderr, "not implemented 245\n"); }
void slaqtr_() { fprintf(stderr, "not implemented 246\n"); }

#endif // DUMMY

} // extern "C"