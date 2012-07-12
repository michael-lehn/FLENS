#ifndef INTEGER
#    ifndef MKL_ILP64
#        define INTEGER int
#    else
#        define INTEGER long
#    endif
#endif

#ifndef FLOAT
#define FLOAT float
#endif

#ifndef DOUBLE
#define DOUBLE double
#endif

#ifndef FLOAT_COMPLEX
#define FLOAT_COMPLEX float 
#endif

#ifndef DOUBLE_COMPLEX
#define DOUBLE_COMPLEX double 
#endif

#ifndef LOGICAL
#define LOGICAL int
#endif

#ifndef UNKNOWN
#define UNKNOWN void
#endif

//-- dbbcsd --------------------------------------------------------------------
void
LAPACK_IMPL(dbbcsd)(const char       *JOBU1,
                    const char       *JOBU2,
                    const char       *JOBV1T,
                    const char       *JOBV2T,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *P,
                    const INTEGER    *Q,
                    DOUBLE           *THETA,
                    DOUBLE           *PHI,
                    DOUBLE           *U1,
                    const INTEGER    *LDU1,
                    DOUBLE           *U2,
                    const INTEGER    *LDU2,
                    DOUBLE           *V1T,
                    const INTEGER    *LDV1T,
                    DOUBLE           *V2T,
                    const INTEGER    *LDV2T,
                    DOUBLE           *B11D,
                    DOUBLE           *B11E,
                    DOUBLE           *B12D,
                    DOUBLE           *B12E,
                    const DOUBLE     *B21D,
                    const DOUBLE     *B21E,
                    const DOUBLE     *B22D,
                    const DOUBLE     *B22E,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dbdsdc --------------------------------------------------------------------
void
LAPACK_IMPL(dbdsdc)(const char       *UPLO,
                    const char       *COMPQ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    DOUBLE           *VT,
                    const INTEGER    *LDVT,
                    DOUBLE           *Q,
                    INTEGER          *IQ,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dbdsqr --------------------------------------------------------------------
void
LAPACK_IMPL(dbdsqr)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NCVT,
                    const INTEGER    *NRU,
                    const INTEGER    *NCC,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *VT,
                    const INTEGER    *LDVT,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- ddisna --------------------------------------------------------------------
void
LAPACK_IMPL(ddisna)(const char       *JOB,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *D,
                    DOUBLE           *SEP,
                    INTEGER          *INFO);

//-- dgbbrd --------------------------------------------------------------------
void
LAPACK_IMPL(dgbbrd)(const char       *VECT,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *NCC,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    DOUBLE           *PT,
                    const INTEGER    *LDPT,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgbcon --------------------------------------------------------------------
void
LAPACK_IMPL(dgbcon)(const char       *NORM,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    const INTEGER    *IPIV,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dgbequ --------------------------------------------------------------------
void
LAPACK_IMPL(dgbequ)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *R,
                    DOUBLE           *C,
                    DOUBLE           *ROWCND,
                    DOUBLE           *COLCND,
                    DOUBLE           *AMAX,
                    INTEGER          *INFO);

//-- dgbequb -------------------------------------------------------------------
void
LAPACK_IMPL(dgbequb)(const INTEGER    *M,
                     const INTEGER    *N,
                     const INTEGER    *KL,
                     const INTEGER    *KU,
                     const DOUBLE     *AB,
                     const INTEGER    *LDAB,
                     DOUBLE           *R,
                     DOUBLE           *C,
                     DOUBLE           *ROWCND,
                     DOUBLE           *COLCND,
                     DOUBLE           *AMAX,
                     INTEGER          *INFO);

//-- dgbrfs --------------------------------------------------------------------
void
LAPACK_IMPL(dgbrfs)(const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    const DOUBLE     *AFB,
                    const INTEGER    *LDAFB,
                    const INTEGER    *IPIV,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dgbsv ---------------------------------------------------------------------
void
LAPACK_IMPL(dgbsv)(const INTEGER        *N,
                   const INTEGER        *KL,
                   const INTEGER        *KU,
                   const INTEGER        *NRHS,
                   DOUBLE               *AB,
                   const INTEGER        *LDAB,
                   INTEGER              *IPIV,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- dgbsvx --------------------------------------------------------------------
void
LAPACK_IMPL(dgbsvx)(const char       *FACT,
                    const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    const INTEGER    *NRHS,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *AFB,
                    const INTEGER    *LDAFB,
                    INTEGER          *IPIV,
                    char             *EQUED,
                    DOUBLE           *R,
                    DOUBLE           *C,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *RCOND,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dgbtf2 --------------------------------------------------------------------
void
LAPACK_IMPL(dgbtf2)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- dgbtrf --------------------------------------------------------------------
void
LAPACK_IMPL(dgbtrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- dgbtrs --------------------------------------------------------------------
void
LAPACK_IMPL(dgbtrs)(const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    const INTEGER    *IPIV,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dgebak --------------------------------------------------------------------
void
LAPACK_IMPL(dgebak)(const char       *JOB,
                    const char       *SIDE,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    const DOUBLE     *SCALE,
                    const INTEGER    *M,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    INTEGER          *INFO);

//-- dgebal --------------------------------------------------------------------
void
LAPACK_IMPL(dgebal)(const char       *JOB,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *ILO,
                    INTEGER          *IHI,
                    DOUBLE           *SCALE,
                    INTEGER          *INFO);

//-- dgebd2 --------------------------------------------------------------------
void
LAPACK_IMPL(dgebd2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *TAUQ,
                    DOUBLE           *TAUP,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgebrd --------------------------------------------------------------------
void
LAPACK_IMPL(dgebrd)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *TAUQ,
                    DOUBLE           *TAUP,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgecon --------------------------------------------------------------------
void
LAPACK_IMPL(dgecon)(const char       *NORM,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dgeequ --------------------------------------------------------------------
void
LAPACK_IMPL(dgeequ)(const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *R,
                    DOUBLE           *C,
                    DOUBLE           *ROWCND,
                    DOUBLE           *COLCND,
                    DOUBLE           *AMAX,
                    INTEGER          *INFO);

//-- dgeequb -------------------------------------------------------------------
void
LAPACK_IMPL(dgeequb)(const INTEGER    *M,
                     const INTEGER    *N,
                     const DOUBLE     *A,
                     const INTEGER    *LDA,
                     DOUBLE           *R,
                     DOUBLE           *C,
                     DOUBLE           *ROWCND,
                     DOUBLE           *COLCND,
                     DOUBLE           *AMAX,
                     INTEGER          *INFO);

//-- dgees ---------------------------------------------------------------------
void
LAPACK_IMPL(dgees)(const char           *JOBVS,
                   const char           *SORT,
                   LOGICAL           (*SELECT)(const DOUBLE *, const DOUBLE *),
                   const INTEGER        *N,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   INTEGER              *SDIM,
                   DOUBLE               *WR,
                   DOUBLE               *WI,
                   DOUBLE               *VS,
                   const INTEGER        *LDVS,
                   DOUBLE               *WORK,
                   const INTEGER        *LWORK,
                   LOGICAL              *BWORK,
                   INTEGER              *INFO);

//-- dgeesx --------------------------------------------------------------------
void
LAPACK_IMPL(dgeesx)(const char       *JOBVS,
                    const char       *SORT,
                    LOGICAL          (*SELECT)(const DOUBLE *, const DOUBLE *),
                    const char       *SENSE,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *SDIM,
                    DOUBLE           *WR,
                    DOUBLE           *WI,
                    DOUBLE           *VS,
                    const INTEGER    *LDVS,
                    DOUBLE           *RCONDE,
                    DOUBLE           *RCONDV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    LOGICAL          *BWORK,
                    INTEGER          *INFO);

//-- dgeev ---------------------------------------------------------------------
void
LAPACK_IMPL(dgeev)(const char           *JOBVL,
                   const char           *JOBVR,
                   const INTEGER        *N,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   DOUBLE               *WR,
                   DOUBLE               *WI,
                   DOUBLE               *VL,
                   const INTEGER        *LDVL,
                   DOUBLE               *VR,
                   const INTEGER        *LDVR,
                   DOUBLE               *WORK,
                   const INTEGER        *LWORK,
                   INTEGER              *INFO);

//-- dgeevx --------------------------------------------------------------------
void
LAPACK_IMPL(dgeevx)(const char       *BALANC,
                    const char       *JOBVL,
                    const char       *JOBVR,
                    const char       *SENSE,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *WR,
                    DOUBLE           *WI,
                    DOUBLE           *VL,
                    const INTEGER    *LDVL,
                    DOUBLE           *VR,
                    const INTEGER    *LDVR,
                    INTEGER          *ILO,
                    INTEGER          *IHI,
                    DOUBLE           *SCALE,
                    DOUBLE           *ABNRM,
                    DOUBLE           *RCONDE,
                    DOUBLE           *RCONDV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dgegs ---------------------------------------------------------------------
void
LAPACK_IMPL(dgegs)(const char           *JOBVSL,
                   const char           *JOBVSR,
                   const INTEGER        *N,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   DOUBLE               *ALPHAR,
                   DOUBLE               *ALPHAI,
                   DOUBLE               *BETA,
                   DOUBLE               *VSL,
                   const INTEGER        *LDVSL,
                   DOUBLE               *VSR,
                   const INTEGER        *LDVSR,
                   DOUBLE               *WORK,
                   const INTEGER        *LWORK,
                   INTEGER              *INFO);

//-- dgegv ---------------------------------------------------------------------
void
LAPACK_IMPL(dgegv)(const char           *JOBVL,
                   const char           *JOBVR,
                   const INTEGER        *N,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   DOUBLE               *ALPHAR,
                   DOUBLE               *ALPHAI,
                   DOUBLE               *BETA,
                   DOUBLE               *VL,
                   const INTEGER        *LDVL,
                   DOUBLE               *VR,
                   const INTEGER        *LDVR,
                   DOUBLE               *WORK,
                   const INTEGER        *LWORK,
                   INTEGER              *INFO);

//-- dgehd2 --------------------------------------------------------------------
void
LAPACK_IMPL(dgehd2)(const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgehrd --------------------------------------------------------------------
void
LAPACK_IMPL(dgehrd)(const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgejsv --------------------------------------------------------------------
void
LAPACK_IMPL(dgejsv)(const char       *JOBA,
                    const char       *JOBU,
                    const char       *JOBV,
                    const char       *JOBR,
                    const char       *JOBT,
                    const char       *JOBP,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *SVA,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dgelq2 --------------------------------------------------------------------
void
LAPACK_IMPL(dgelq2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgelqf --------------------------------------------------------------------
void
LAPACK_IMPL(dgelqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgels ---------------------------------------------------------------------
void
LAPACK_IMPL(dgels)(const char           *TRANS,
                   const INTEGER        *M,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   DOUBLE               *WORK,
                   const INTEGER        *LWORK,
                   INTEGER              *INFO);

//-- dgelsd --------------------------------------------------------------------
void
LAPACK_IMPL(dgelsd)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *S,
                    const DOUBLE     *RCOND,
                    INTEGER          *RANK,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dgelss --------------------------------------------------------------------
void
LAPACK_IMPL(dgelss)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *S,
                    const DOUBLE     *RCOND,
                    INTEGER          *RANK,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgelsx --------------------------------------------------------------------
void
LAPACK_IMPL(dgelsx)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *JPVT,
                    const DOUBLE     *RCOND,
                    INTEGER          *RANK,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgelsy --------------------------------------------------------------------
void
LAPACK_IMPL(dgelsy)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *JPVT,
                    const DOUBLE     *RCOND,
                    INTEGER          *RANK,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgeql2 --------------------------------------------------------------------
void
LAPACK_IMPL(dgeql2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgeqlf --------------------------------------------------------------------
void
LAPACK_IMPL(dgeqlf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgeqp3 --------------------------------------------------------------------
void
LAPACK_IMPL(dgeqp3)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *JPVT,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgeqpf --------------------------------------------------------------------
void
LAPACK_IMPL(dgeqpf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *JPVT,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgeqr2 --------------------------------------------------------------------
void
LAPACK_IMPL(dgeqr2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgeqr2p -------------------------------------------------------------------
void
LAPACK_IMPL(dgeqr2p)(const INTEGER    *M,
                     const INTEGER    *N,
                     DOUBLE           *A,
                     const INTEGER    *LDA,
                     DOUBLE           *TAU,
                     DOUBLE           *WORK,
                     INTEGER          *INFO);

//-- dgeqrf --------------------------------------------------------------------
void
LAPACK_IMPL(dgeqrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgeqrfp -------------------------------------------------------------------
void
LAPACK_IMPL(dgeqrfp)(const INTEGER    *M,
                     const INTEGER    *N,
                     DOUBLE           *A,
                     const INTEGER    *LDA,
                     DOUBLE           *TAU,
                     DOUBLE           *WORK,
                     const INTEGER    *LWORK,
                     INTEGER          *INFO);

//-- dgerfs --------------------------------------------------------------------
void
LAPACK_IMPL(dgerfs)(const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *AF,
                    const INTEGER    *LDAF,
                    const INTEGER    *IPIV,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dgerq2 --------------------------------------------------------------------
void
LAPACK_IMPL(dgerq2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgerqf --------------------------------------------------------------------
void
LAPACK_IMPL(dgerqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgesc2 --------------------------------------------------------------------
void
LAPACK_IMPL(dgesc2)(const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *RHS,
                    const INTEGER    *IPIV,
                    const INTEGER    *JPIV,
                    DOUBLE           *SCALE);

//-- dgesdd --------------------------------------------------------------------
void
LAPACK_IMPL(dgesdd)(const char       *JOBZ,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *S,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    DOUBLE           *VT,
                    const INTEGER    *LDVT,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dgesv ---------------------------------------------------------------------
void
LAPACK_IMPL(dgesv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   INTEGER              *IPIV,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- dgesvd --------------------------------------------------------------------
void
LAPACK_IMPL(dgesvd)(const char       *JOBU,
                    const char       *JOBVT,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *S,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    DOUBLE           *VT,
                    const INTEGER    *LDVT,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgesvj --------------------------------------------------------------------
void
LAPACK_IMPL(dgesvj)(const char       *JOBA,
                    const char       *JOBU,
                    const char       *JOBV,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *SVA,
                    const INTEGER    *MV,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgesvx --------------------------------------------------------------------
void
LAPACK_IMPL(dgesvx)(const char       *FACT,
                    const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *AF,
                    const INTEGER    *LDAF,
                    INTEGER          *IPIV,
                    char             *EQUED,
                    DOUBLE           *R,
                    DOUBLE           *C,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *RCOND,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dgetc2 --------------------------------------------------------------------
void
LAPACK_IMPL(dgetc2)(const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *JPIV,
                    INTEGER          *INFO);

//-- dgetf2 --------------------------------------------------------------------
void
LAPACK_IMPL(dgetf2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- dgetrf --------------------------------------------------------------------
void
LAPACK_IMPL(dgetrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- dgetri --------------------------------------------------------------------
void
LAPACK_IMPL(dgetri)(const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgetrs --------------------------------------------------------------------
void
LAPACK_IMPL(dgetrs)(const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dggbak --------------------------------------------------------------------
void
LAPACK_IMPL(dggbak)(const char       *JOB,
                    const char       *SIDE,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    const DOUBLE     *LSCALE,
                    const DOUBLE     *RSCALE,
                    const INTEGER    *M,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    INTEGER          *INFO);

//-- dggbal --------------------------------------------------------------------
void
LAPACK_IMPL(dggbal)(const char       *JOB,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *ILO,
                    INTEGER          *IHI,
                    DOUBLE           *LSCALE,
                    DOUBLE           *RSCALE,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgges ---------------------------------------------------------------------
void
LAPACK_IMPL(dgges)(const char           *JOBVSL,
                   const char           *JOBVSR,
                   const char           *SORT,
                   const LOGICAL        *SELCTG,
                   const INTEGER        *N,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *SDIM,
                   DOUBLE               *ALPHAR,
                   DOUBLE               *ALPHAI,
                   DOUBLE               *BETA,
                   DOUBLE               *VSL,
                   const INTEGER        *LDVSL,
                   DOUBLE               *VSR,
                   const INTEGER        *LDVSR,
                   DOUBLE               *WORK,
                   const INTEGER        *LWORK,
                   LOGICAL              *BWORK,
                   INTEGER              *INFO);

//-- dggesx --------------------------------------------------------------------
void
LAPACK_IMPL(dggesx)(const char       *JOBVSL,
                    const char       *JOBVSR,
                    const char       *SORT,
                    const LOGICAL    *SELCTG,
                    const char       *SENSE,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *SDIM,
                    DOUBLE           *ALPHAR,
                    DOUBLE           *ALPHAI,
                    DOUBLE           *BETA,
                    DOUBLE           *VSL,
                    const INTEGER    *LDVSL,
                    DOUBLE           *VSR,
                    const INTEGER    *LDVSR,
                    DOUBLE           *RCONDE,
                    DOUBLE           *RCONDV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    LOGICAL          *BWORK,
                    INTEGER          *INFO);

//-- dggev ---------------------------------------------------------------------
void
LAPACK_IMPL(dggev)(const char           *JOBVL,
                   const char           *JOBVR,
                   const INTEGER        *N,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   DOUBLE               *ALPHAR,
                   DOUBLE               *ALPHAI,
                   DOUBLE               *BETA,
                   DOUBLE               *VL,
                   const INTEGER        *LDVL,
                   DOUBLE               *VR,
                   const INTEGER        *LDVR,
                   DOUBLE               *WORK,
                   const INTEGER        *LWORK,
                   INTEGER              *INFO);

//-- dggevx --------------------------------------------------------------------
void
LAPACK_IMPL(dggevx)(const char       *BALANC,
                    const char       *JOBVL,
                    const char       *JOBVR,
                    const char       *SENSE,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *ALPHAR,
                    DOUBLE           *ALPHAI,
                    DOUBLE           *BETA,
                    DOUBLE           *VL,
                    const INTEGER    *LDVL,
                    DOUBLE           *VR,
                    const INTEGER    *LDVR,
                    INTEGER          *ILO,
                    INTEGER          *IHI,
                    DOUBLE           *LSCALE,
                    DOUBLE           *RSCALE,
                    DOUBLE           *ABNRM,
                    DOUBLE           *BBNRM,
                    DOUBLE           *RCONDE,
                    DOUBLE           *RCONDV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    LOGICAL          *BWORK,
                    INTEGER          *INFO);

//-- dggglm --------------------------------------------------------------------
void
LAPACK_IMPL(dggglm)(const INTEGER    *N,
                    const INTEGER    *M,
                    const INTEGER    *P,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *D,
                    DOUBLE           *X,
                    DOUBLE           *Y,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgghrd --------------------------------------------------------------------
void
LAPACK_IMPL(dgghrd)(const char       *COMPQ,
                    const char       *COMPZ,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *INFO);

//-- dgglse --------------------------------------------------------------------
void
LAPACK_IMPL(dgglse)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *P,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *C,
                    DOUBLE           *D,
                    DOUBLE           *X,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dggqrf --------------------------------------------------------------------
void
LAPACK_IMPL(dggqrf)(const INTEGER    *N,
                    const INTEGER    *M,
                    const INTEGER    *P,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAUA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *TAUB,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dggrqf --------------------------------------------------------------------
void
LAPACK_IMPL(dggrqf)(const INTEGER    *M,
                    const INTEGER    *P,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAUA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *TAUB,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dggsvd --------------------------------------------------------------------
void
LAPACK_IMPL(dggsvd)(const char       *JOBU,
                    const char       *JOBV,
                    const char       *JOBQ,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *P,
                    INTEGER          *K,
                    INTEGER          *L,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *ALPHA,
                    DOUBLE           *BETA,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dggsvp --------------------------------------------------------------------
void
LAPACK_IMPL(dggsvp)(const char       *JOBU,
                    const char       *JOBV,
                    const char       *JOBQ,
                    const INTEGER    *M,
                    const INTEGER    *P,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    const DOUBLE     *TOLA,
                    const DOUBLE     *TOLB,
                    INTEGER          *K,
                    INTEGER          *L,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    INTEGER          *IWORK,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgsvj0 --------------------------------------------------------------------
void
LAPACK_IMPL(dgsvj0)(const char       *JOBV,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *SVA,
                    const INTEGER    *MV,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    const DOUBLE     *EPS,
                    const DOUBLE     *SFMIN,
                    const DOUBLE     *TOL,
                    const INTEGER    *NSWEEP,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgsvj1 --------------------------------------------------------------------
void
LAPACK_IMPL(dgsvj1)(const char       *JOBV,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *N1,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *SVA,
                    const INTEGER    *MV,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    const DOUBLE     *EPS,
                    const DOUBLE     *SFMIN,
                    const DOUBLE     *TOL,
                    const INTEGER    *NSWEEP,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgtcon --------------------------------------------------------------------
void
LAPACK_IMPL(dgtcon)(const char       *NORM,
                    const INTEGER    *N,
                    const DOUBLE     *DL,
                    const DOUBLE     *D,
                    const DOUBLE     *DU,
                    const DOUBLE     *DU2,
                    const INTEGER    *IPIV,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dgtrfs --------------------------------------------------------------------
void
LAPACK_IMPL(dgtrfs)(const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *DL,
                    const DOUBLE     *D,
                    const DOUBLE     *DU,
                    const DOUBLE     *DLF,
                    const DOUBLE     *DF,
                    const DOUBLE     *DUF,
                    const DOUBLE     *DU2,
                    const INTEGER    *IPIV,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dgtsv ---------------------------------------------------------------------
void
LAPACK_IMPL(dgtsv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *DL,
                   DOUBLE               *D,
                   DOUBLE               *DU,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- dgtsvx --------------------------------------------------------------------
void
LAPACK_IMPL(dgtsvx)(const char       *FACT,
                    const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *DL,
                    const DOUBLE     *D,
                    const DOUBLE     *DU,
                    DOUBLE           *DLF,
                    DOUBLE           *DF,
                    DOUBLE           *DUF,
                    DOUBLE           *DU2,
                    INTEGER          *IPIV,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *RCOND,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dgttrf --------------------------------------------------------------------
void
LAPACK_IMPL(dgttrf)(const INTEGER    *N,
                    DOUBLE           *DL,
                    DOUBLE           *D,
                    DOUBLE           *DU,
                    DOUBLE           *DU2,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- dgttrs --------------------------------------------------------------------
void
LAPACK_IMPL(dgttrs)(const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *DL,
                    const DOUBLE     *D,
                    const DOUBLE     *DU,
                    const DOUBLE     *DU2,
                    const INTEGER    *IPIV,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dgtts2 --------------------------------------------------------------------
void
LAPACK_IMPL(dgtts2)(const INTEGER    *ITRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *DL,
                    const DOUBLE     *D,
                    const DOUBLE     *DU,
                    const DOUBLE     *DU2,
                    const INTEGER    *IPIV,
                    DOUBLE           *B,
                    const INTEGER    *LDB);

//-- dhgeqz --------------------------------------------------------------------
void
LAPACK_IMPL(dhgeqz)(const char       *JOB,
                    const char       *COMPQ,
                    const char       *COMPZ,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE           *H,
                    const INTEGER    *LDH,
                    DOUBLE           *T,
                    const INTEGER    *LDT,
                    DOUBLE           *ALPHAR,
                    DOUBLE           *ALPHAI,
                    DOUBLE           *BETA,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dhsein --------------------------------------------------------------------
void
LAPACK_IMPL(dhsein)(const char       *SIDE,
                    const char       *EIGSRC,
                    const char       *INITV,
                    LOGICAL          *SELECT,
                    const INTEGER    *N,
                    const DOUBLE     *H,
                    const INTEGER    *LDH,
                    DOUBLE           *WR,
                    const DOUBLE     *WI,
                    DOUBLE           *VL,
                    const INTEGER    *LDVL,
                    DOUBLE           *VR,
                    const INTEGER    *LDVR,
                    const INTEGER    *MM,
                    INTEGER          *M,
                    DOUBLE           *WORK,
                    INTEGER          *IFAILL,
                    INTEGER          *IFAILR,
                    INTEGER          *INFO);

//-- dhseqr --------------------------------------------------------------------
void
LAPACK_IMPL(dhseqr)(const char       *JOB,
                    const char       *COMPZ,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE           *H,
                    const INTEGER    *LDH,
                    DOUBLE           *WR,
                    DOUBLE           *WI,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- disnan --------------------------------------------------------------------
LOGICAL
LAPACK_IMPL(disnan)(const DOUBLE     *DIN);

//-- dla_gbamv -----------------------------------------------------------------
void
LAPACK_IMPL(dla_gbamv)(const INTEGER        *TRANS,
                       const INTEGER        *M,
                       const INTEGER        *N,
                       const INTEGER        *KL,
                       const INTEGER        *KU,
                       const DOUBLE         *ALPHA,
                       const DOUBLE         *AB,
                       const INTEGER        *LDAB,
                       const DOUBLE         *X,
                       const INTEGER        *INCX,
                       const DOUBLE         *BETA,
                       DOUBLE               *Y,
                       const INTEGER        *INCY);

//-- dla_gbrcond ---------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dla_gbrcond)(const char       *TRANS,
                         const INTEGER    *N,
                         const INTEGER    *KL,
                         const INTEGER    *KU,
                         const DOUBLE     *AB,
                         const INTEGER    *LDAB,
                         const DOUBLE     *AFB,
                         const INTEGER    *LDAFB,
                         const INTEGER    *IPIV,
                         const INTEGER    *CMODE,
                         const DOUBLE     *C,
                         INTEGER          *INFO,
                         const DOUBLE     *WORK,
                         const INTEGER    *IWORK);

//-- dla_gbrpvgrw --------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dla_gbrpvgrw)(const INTEGER    *N,
                          const INTEGER    *KL,
                          const INTEGER    *KU,
                          const INTEGER    *NCOLS,
                          const DOUBLE     *AB,
                          const INTEGER    *LDAB,
                          const DOUBLE     *AFB,
                          const INTEGER    *LDAFB);

//-- dla_geamv -----------------------------------------------------------------
void
LAPACK_IMPL(dla_geamv)(const INTEGER        *TRANS,
                       const INTEGER        *M,
                       const INTEGER        *N,
                       const DOUBLE         *ALPHA,
                       const DOUBLE         *A,
                       const INTEGER        *LDA,
                       const DOUBLE         *X,
                       const INTEGER        *INCX,
                       const DOUBLE         *BETA,
                       DOUBLE               *Y,
                       const INTEGER        *INCY);

//-- dla_gercond ---------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dla_gercond)(const char       *TRANS,
                         const INTEGER    *N,
                         const DOUBLE     *A,
                         const INTEGER    *LDA,
                         const DOUBLE     *AF,
                         const INTEGER    *LDAF,
                         const INTEGER    *IPIV,
                         const INTEGER    *CMODE,
                         const DOUBLE     *C,
                         INTEGER          *INFO,
                         const DOUBLE     *WORK,
                         const INTEGER    *IWORK);

//-- dla_lin_berr --------------------------------------------------------------
void
LAPACK_IMPL(dla_lin_berr)(const INTEGER    *N,
                          const INTEGER    *NZ,
                          const INTEGER    *NRHS,
                          const DOUBLE     *RES,
                          const DOUBLE     *AYB,
                          DOUBLE           *BERR);

//-- dla_porcond ---------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dla_porcond)(const char       *UPLO,
                         const INTEGER    *N,
                         const DOUBLE     *A,
                         const INTEGER    *LDA,
                         const DOUBLE     *AF,
                         const INTEGER    *LDAF,
                         const INTEGER    *CMODE,
                         const DOUBLE     *C,
                         INTEGER          *INFO,
                         const DOUBLE     *WORK,
                         const INTEGER    *IWORK);

//-- dla_porpvgrw --------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dla_porpvgrw)(const char       *UPLO,
                          const INTEGER    *NCOLS,
                          const DOUBLE     *A,
                          const INTEGER    *LDA,
                          const DOUBLE     *AF,
                          const INTEGER    *LDAF,
                          const DOUBLE     *WORK);

//-- dla_rpvgrw ----------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dla_rpvgrw)(const INTEGER    *N,
                        const INTEGER    *NCOLS,
                        const DOUBLE     *A,
                        const INTEGER    *LDA,
                        const DOUBLE     *AF,
                        const INTEGER    *LDAF);

//-- dla_syamv -----------------------------------------------------------------
void
LAPACK_IMPL(dla_syamv)(const INTEGER        *UPLO,
                       const INTEGER        *N,
                       const DOUBLE         *ALPHA,
                       const DOUBLE         *A,
                       const INTEGER        *LDA,
                       const DOUBLE         *X,
                       const INTEGER        *INCX,
                       const DOUBLE         *BETA,
                       DOUBLE               *Y,
                       const INTEGER        *INCY);

//-- dla_syrcond ---------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dla_syrcond)(const char       *UPLO,
                         const INTEGER    *N,
                         const DOUBLE     *A,
                         const INTEGER    *LDA,
                         const DOUBLE     *AF,
                         const INTEGER    *LDAF,
                         const INTEGER    *IPIV,
                         const INTEGER    *CMODE,
                         const DOUBLE     *C,
                         INTEGER          *INFO,
                         const DOUBLE     *WORK,
                         const INTEGER    *IWORK);

//-- dla_syrpvgrw --------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dla_syrpvgrw)(const char       *UPLO,
                          const INTEGER    *N,
                          const INTEGER    *INFO,
                          const DOUBLE     *A,
                          const INTEGER    *LDA,
                          const DOUBLE     *AF,
                          const INTEGER    *LDAF,
                          const INTEGER    *IPIV,
                          const DOUBLE     *WORK);

//-- dla_wwaddw ----------------------------------------------------------------
void
LAPACK_IMPL(dla_wwaddw)(const INTEGER    *N,
                        DOUBLE           *X,
                        DOUBLE           *Y,
                        const DOUBLE     *W);

//-- dlabad --------------------------------------------------------------------
void
LAPACK_IMPL(dlabad)(DOUBLE   *SMALL,
                    DOUBLE   *LARGE);

//-- dlabrd --------------------------------------------------------------------
void
LAPACK_IMPL(dlabrd)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *NB,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *TAUQ,
                    DOUBLE           *TAUP,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *Y,
                    const INTEGER    *LDY);

//-- dlacn2 --------------------------------------------------------------------
void
LAPACK_IMPL(dlacn2)(const INTEGER    *N,
                    DOUBLE           *V,
                    DOUBLE           *X,
                    INTEGER          *ISGN,
                    DOUBLE           *EST,
                    INTEGER          *KASE,
                    INTEGER          *ISAVE);

//-- dlacon --------------------------------------------------------------------
void
LAPACK_IMPL(dlacon)(const INTEGER    *N,
                    DOUBLE           *V,
                    DOUBLE           *X,
                    INTEGER          *ISGN,
                    DOUBLE           *EST,
                    INTEGER          *KASE);

//-- dlacpy --------------------------------------------------------------------
void
LAPACK_IMPL(dlacpy)(const char       *UPLO,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB);

//-- dladiv --------------------------------------------------------------------
void
LAPACK_IMPL(dladiv)(const DOUBLE     *A,
                    const DOUBLE     *B,
                    const DOUBLE     *C,
                    const DOUBLE     *D,
                    DOUBLE           *P,
                    DOUBLE           *Q);

//-- dlae2 ---------------------------------------------------------------------
void
LAPACK_IMPL(dlae2)(const DOUBLE     *A,
                   const DOUBLE     *B,
                   const DOUBLE     *C,
                   DOUBLE           *RT1,
                   DOUBLE           *RT2);

//-- dlaebz --------------------------------------------------------------------
void
LAPACK_IMPL(dlaebz)(const INTEGER    *IJOB,
                    const INTEGER    *NITMAX,
                    const INTEGER    *N,
                    const INTEGER    *MMAX,
                    const INTEGER    *MINP,
                    const INTEGER    *NBMIN,
                    const DOUBLE     *ABSTOL,
                    const DOUBLE     *RELTOL,
                    const DOUBLE     *PIVMIN,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    const DOUBLE     *E2,
                    INTEGER          *NVAL,
                    DOUBLE           *AB,
                    DOUBLE           *C,
                    INTEGER          *MOUT,
                    INTEGER          *NAB,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dlaed0 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaed0)(const INTEGER    *ICOMPQ,
                    const INTEGER    *QSIZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    const DOUBLE     *E,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    DOUBLE           *QSTORE,
                    const INTEGER    *LDQS,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dlaed1 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaed1)(const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    INTEGER          *INDXQ,
                    const DOUBLE     *RHO,
                    const INTEGER    *CUTPNT,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dlaed2 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaed2)(INTEGER          *K,
                    const INTEGER    *N,
                    const INTEGER    *N1,
                    DOUBLE           *D,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    INTEGER          *INDXQ,
                    DOUBLE           *RHO,
                    const DOUBLE     *Z,
                    DOUBLE           *DLAMDA,
                    DOUBLE           *W,
                    DOUBLE           *Q2,
                    INTEGER          *INDX,
                    INTEGER          *INDXC,
                    INTEGER          *INDXP,
                    INTEGER          *COLTYP,
                    INTEGER          *INFO);

//-- dlaed3 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaed3)(const INTEGER    *K,
                    const INTEGER    *N,
                    const INTEGER    *N1,
                    DOUBLE           *D,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    const DOUBLE     *RHO,
                    DOUBLE           *DLAMDA,
                    const DOUBLE     *Q2,
                    const INTEGER    *INDX,
                    const INTEGER    *CTOT,
                    DOUBLE           *W,
                    DOUBLE           *S,
                    INTEGER          *INFO);

//-- dlaed4 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaed4)(const INTEGER    *N,
                    const INTEGER    *I,
                    const DOUBLE     *D,
                    const DOUBLE     *Z,
                    DOUBLE           *DELTA,
                    const DOUBLE     *RHO,
                    DOUBLE           *DLAM,
                    INTEGER          *INFO);

//-- dlaed5 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaed5)(const INTEGER    *I,
                    const DOUBLE     *D,
                    const DOUBLE     *Z,
                    DOUBLE           *DELTA,
                    const DOUBLE     *RHO,
                    DOUBLE           *DLAM);

//-- dlaed6 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaed6)(const INTEGER    *KNITER,
                    const LOGICAL    *ORGATI,
                    const DOUBLE     *RHO,
                    const DOUBLE     *D,
                    const DOUBLE     *Z,
                    const DOUBLE     *FINIT,
                    DOUBLE           *TAU,
                    INTEGER          *INFO);

//-- dlaed7 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaed7)(const INTEGER    *ICOMPQ,
                    const INTEGER    *N,
                    const INTEGER    *QSIZ,
                    const INTEGER    *TLVLS,
                    const INTEGER    *CURLVL,
                    const INTEGER    *CURPBM,
                    DOUBLE           *D,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    INTEGER          *INDXQ,
                    const DOUBLE     *RHO,
                    const INTEGER    *CUTPNT,
                    DOUBLE           *QSTORE,
                    INTEGER          *QPTR,
                    const INTEGER    *PRMPTR,
                    const INTEGER    *PERM,
                    const INTEGER    *GIVPTR,
                    const INTEGER    *GIVCOL,
                    const DOUBLE     *GIVNUM,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dlaed8 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaed8)(const INTEGER    *ICOMPQ,
                    INTEGER          *K,
                    const INTEGER    *N,
                    const INTEGER    *QSIZ,
                    DOUBLE           *D,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    const INTEGER    *INDXQ,
                    DOUBLE           *RHO,
                    const INTEGER    *CUTPNT,
                    const DOUBLE     *Z,
                    DOUBLE           *DLAMDA,
                    DOUBLE           *Q2,
                    const INTEGER    *LDQ2,
                    DOUBLE           *W,
                    INTEGER          *PERM,
                    INTEGER          *GIVPTR,
                    INTEGER          *GIVCOL,
                    DOUBLE           *GIVNUM,
                    INTEGER          *INDXP,
                    INTEGER          *INDX,
                    INTEGER          *INFO);

//-- dlaed9 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaed9)(const INTEGER    *K,
                    const INTEGER    *KSTART,
                    const INTEGER    *KSTOP,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    const DOUBLE     *RHO,
                    const DOUBLE     *DLAMDA,
                    const DOUBLE     *W,
                    DOUBLE           *S,
                    const INTEGER    *LDS,
                    INTEGER          *INFO);

//-- dlaeda --------------------------------------------------------------------
void
LAPACK_IMPL(dlaeda)(const INTEGER    *N,
                    const INTEGER    *TLVLS,
                    const INTEGER    *CURLVL,
                    const INTEGER    *CURPBM,
                    const INTEGER    *PRMPTR,
                    const INTEGER    *PERM,
                    const INTEGER    *GIVPTR,
                    const INTEGER    *GIVCOL,
                    const DOUBLE     *GIVNUM,
                    const DOUBLE     *Q,
                    const INTEGER    *QPTR,
                    DOUBLE           *Z,
                    DOUBLE           *ZTEMP,
                    INTEGER          *INFO);

//-- dlaein --------------------------------------------------------------------
void
LAPACK_IMPL(dlaein)(const LOGICAL    *RIGHTV,
                    const LOGICAL    *NOINIT,
                    const INTEGER    *N,
                    const DOUBLE     *H,
                    const INTEGER    *LDH,
                    const DOUBLE     *WR,
                    const DOUBLE     *WI,
                    DOUBLE           *VR,
                    DOUBLE           *VI,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *WORK,
                    const DOUBLE     *EPS3,
                    const DOUBLE     *SMLNUM,
                    const DOUBLE     *BIGNUM,
                    INTEGER          *INFO);

//-- dlaev2 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaev2)(const DOUBLE     *A,
                    const DOUBLE     *B,
                    const DOUBLE     *C,
                    DOUBLE           *RT1,
                    DOUBLE           *RT2,
                    DOUBLE           *CS1,
                    DOUBLE           *SN1);

//-- dlaexc --------------------------------------------------------------------
void
LAPACK_IMPL(dlaexc)(const LOGICAL    *WANTQ,
                    const INTEGER    *N,
                    DOUBLE           *T,
                    const INTEGER    *LDT,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    const INTEGER    *J1,
                    const INTEGER    *N1,
                    const INTEGER    *N2,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dlag2 ---------------------------------------------------------------------
void
LAPACK_IMPL(dlag2)(const DOUBLE         *A,
                   const INTEGER        *LDA,
                   const DOUBLE         *B,
                   const INTEGER        *LDB,
                   const DOUBLE         *SAFMIN,
                   DOUBLE               *SCALE1,
                   DOUBLE               *SCALE2,
                   DOUBLE               *WR1,
                   DOUBLE               *WR2,
                   DOUBLE               *WI);

//-- dlag2s --------------------------------------------------------------------
void
LAPACK_IMPL(dlag2s)(const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    FLOAT            *SA,
                    const INTEGER    *LDSA,
                    INTEGER          *INFO);

//-- dlags2 --------------------------------------------------------------------
void
LAPACK_IMPL(dlags2)(const LOGICAL    *UPPER,
                    const DOUBLE     *A1,
                    const DOUBLE     *A2,
                    const DOUBLE     *A3,
                    const DOUBLE     *B1,
                    const DOUBLE     *B2,
                    const DOUBLE     *B3,
                    DOUBLE           *CSU,
                    DOUBLE           *SNU,
                    DOUBLE           *CSV,
                    DOUBLE           *SNV,
                    DOUBLE           *CSQ,
                    DOUBLE           *SNQ);

//-- dlagtf --------------------------------------------------------------------
void
LAPACK_IMPL(dlagtf)(const INTEGER    *N,
                    DOUBLE           *A,
                    const DOUBLE     *LAMBDA,
                    DOUBLE           *B,
                    DOUBLE           *C,
                    const DOUBLE     *TOL,
                    DOUBLE           *D,
                    INTEGER          *IN,
                    INTEGER          *INFO);

//-- dlagtm --------------------------------------------------------------------
void
LAPACK_IMPL(dlagtm)(const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *ALPHA,
                    const DOUBLE     *DL,
                    const DOUBLE     *D,
                    const DOUBLE     *DU,
                    const DOUBLE     *X,
                    const INTEGER    *LDX,
                    const DOUBLE     *BETA,
                    DOUBLE           *B,
                    const INTEGER    *LDB);

//-- dlagts --------------------------------------------------------------------
void
LAPACK_IMPL(dlagts)(const INTEGER    *JOB,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const DOUBLE     *B,
                    const DOUBLE     *C,
                    const DOUBLE     *D,
                    const INTEGER    *IN,
                    DOUBLE           *Y,
                    DOUBLE           *TOL,
                    INTEGER          *INFO);

//-- dlagv2 --------------------------------------------------------------------
void
LAPACK_IMPL(dlagv2)(DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *ALPHAR,
                    DOUBLE           *ALPHAI,
                    DOUBLE           *BETA,
                    DOUBLE           *CSL,
                    DOUBLE           *SNL,
                    DOUBLE           *CSR,
                    DOUBLE           *SNR);

//-- dlahqr --------------------------------------------------------------------
void
LAPACK_IMPL(dlahqr)(const LOGICAL    *WANTT,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE           *H,
                    const INTEGER    *LDH,
                    DOUBLE           *WR,
                    DOUBLE           *WI,
                    const INTEGER    *ILOZ,
                    const INTEGER    *IHIZ,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *INFO);

//-- dlahr2 --------------------------------------------------------------------
void
LAPACK_IMPL(dlahr2)(const INTEGER    *N,
                    const INTEGER    *K,
                    const INTEGER    *NB,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *T,
                    const INTEGER    *LDT,
                    DOUBLE           *Y,
                    const INTEGER    *LDY);

//-- dlahrd --------------------------------------------------------------------
void
LAPACK_IMPL(dlahrd)(const INTEGER    *N,
                    const INTEGER    *K,
                    const INTEGER    *NB,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *T,
                    const INTEGER    *LDT,
                    DOUBLE           *Y,
                    const INTEGER    *LDY);

//-- dlaic1 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaic1)(const INTEGER    *JOB,
                    const INTEGER    *J,
                    const DOUBLE     *X,
                    const DOUBLE     *SEST,
                    const DOUBLE     *W,
                    const DOUBLE     *GAMMA,
                    DOUBLE           *SESTPR,
                    DOUBLE           *S,
                    DOUBLE           *C);

//-- dlaisnan ------------------------------------------------------------------
LOGICAL
LAPACK_IMPL(dlaisnan)(const DOUBLE     *DIN1,
                      const DOUBLE     *DIN2);

//-- dlaln2 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaln2)(const LOGICAL    *LTRANS,
                    const INTEGER    *NA,
                    const INTEGER    *NW,
                    const DOUBLE     *SMIN,
                    const DOUBLE     *CA,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *D1,
                    const DOUBLE     *D2,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    const DOUBLE     *WR,
                    const DOUBLE     *WI,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *SCALE,
                    DOUBLE           *XNORM,
                    INTEGER          *INFO);

//-- dlals0 --------------------------------------------------------------------
void
LAPACK_IMPL(dlals0)(const INTEGER    *ICOMPQ,
                    const INTEGER    *NL,
                    const INTEGER    *NR,
                    const INTEGER    *SQRE,
                    const INTEGER    *NRHS,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *BX,
                    const INTEGER    *LDBX,
                    const INTEGER    *PERM,
                    const INTEGER    *GIVPTR,
                    const INTEGER    *GIVCOL,
                    const INTEGER    *LDGCOL,
                    const DOUBLE     *GIVNUM,
                    const INTEGER    *LDGNUM,
                    const DOUBLE     *POLES,
                    const DOUBLE     *DIFL,
                    const DOUBLE     *DIFR,
                    const DOUBLE     *Z,
                    const INTEGER    *K,
                    const DOUBLE     *C,
                    const DOUBLE     *S,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dlalsa --------------------------------------------------------------------
void
LAPACK_IMPL(dlalsa)(const INTEGER    *ICOMPQ,
                    const INTEGER    *SMLSIZ,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *BX,
                    const INTEGER    *LDBX,
                    const DOUBLE     *U,
                    const INTEGER    *LDU,
                    const DOUBLE     *VT,
                    const INTEGER    *K,
                    const DOUBLE     *DIFL,
                    const DOUBLE     *DIFR,
                    const DOUBLE     *Z,
                    const DOUBLE     *POLES,
                    const INTEGER    *GIVPTR,
                    const INTEGER    *GIVCOL,
                    const INTEGER    *LDGCOL,
                    const INTEGER    *PERM,
                    const DOUBLE     *GIVNUM,
                    const DOUBLE     *C,
                    const DOUBLE     *S,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dlalsd --------------------------------------------------------------------
void
LAPACK_IMPL(dlalsd)(const char       *UPLO,
                    const INTEGER    *SMLSIZ,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    const DOUBLE     *RCOND,
                    INTEGER          *RANK,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dlamch --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dlamch)(const char   *CMACH);

//-- dlamrg --------------------------------------------------------------------
void
LAPACK_IMPL(dlamrg)(const INTEGER    *N1,
                    const INTEGER    *N2,
                    const DOUBLE     *A,
                    const INTEGER    *DTRD1,
                    const INTEGER    *DTRD2,
                    INTEGER          *INDEX);

//-- dlaneg --------------------------------------------------------------------
INTEGER
LAPACK_IMPL(dlaneg)(const INTEGER    *N,
                    const DOUBLE     *D,
                    const DOUBLE     *LLD,
                    const DOUBLE     *SIGMA,
                    const DOUBLE     *PIVMIN,
                    const INTEGER    *R);

//-- dlangb --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dlangb)(const char       *NORM,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *WORK);

//-- dlange --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dlange)(const char       *NORM,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *WORK);

//-- dlangt --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dlangt)(const char       *NORM,
                    const INTEGER    *N,
                    const DOUBLE     *DL,
                    const DOUBLE     *D,
                    const DOUBLE     *DU);

//-- dlanhs --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dlanhs)(const char       *NORM,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *WORK);

//-- dlansb --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dlansb)(const char       *NORM,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *WORK);

//-- dlansf --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dlansf)(const char       *NORM,
                    const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    DOUBLE           *WORK);

//-- dlansp --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dlansp)(const char       *NORM,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *WORK);

//-- dlanst --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dlanst)(const char       *NORM,
                    const INTEGER    *N,
                    const DOUBLE     *D,
                    const DOUBLE     *E);

//-- dlansy --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dlansy)(const char       *NORM,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *WORK);

//-- dlantb --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dlantb)(const char       *NORM,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *WORK);

//-- dlantp --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dlantp)(const char       *NORM,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *WORK);

//-- dlantr --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dlantr)(const char       *NORM,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *WORK);

//-- dlanv2 --------------------------------------------------------------------
void
LAPACK_IMPL(dlanv2)(DOUBLE   *A,
                    DOUBLE   *B,
                    DOUBLE   *C,
                    DOUBLE   *D,
                    DOUBLE   *RT1R,
                    DOUBLE   *RT1I,
                    DOUBLE   *RT2R,
                    DOUBLE   *RT2I,
                    DOUBLE   *CS,
                    DOUBLE   *SN);

//-- dlapll --------------------------------------------------------------------
void
LAPACK_IMPL(dlapll)(const INTEGER    *N,
                    DOUBLE           *X,
                    const INTEGER    *INCX,
                    DOUBLE           *Y,
                    const INTEGER    *INCY,
                    DOUBLE           *SSMIN);

//-- dlapmr --------------------------------------------------------------------
void
LAPACK_IMPL(dlapmr)(const LOGICAL    *FORWRD,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    INTEGER          *K);

//-- dlapmt --------------------------------------------------------------------
void
LAPACK_IMPL(dlapmt)(const LOGICAL    *FORWRD,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    INTEGER          *K);

//-- dlapy2 --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dlapy2)(const DOUBLE     *X,
                    const DOUBLE     *Y);

//-- dlapy3 --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dlapy3)(const DOUBLE     *X,
                    const DOUBLE     *Y,
                    const DOUBLE     *Z);

//-- dlaqgb --------------------------------------------------------------------
void
LAPACK_IMPL(dlaqgb)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    const DOUBLE     *R,
                    const DOUBLE     *C,
                    const DOUBLE     *ROWCND,
                    const DOUBLE     *COLCND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- dlaqge --------------------------------------------------------------------
void
LAPACK_IMPL(dlaqge)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *R,
                    const DOUBLE     *C,
                    const DOUBLE     *ROWCND,
                    const DOUBLE     *COLCND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- dlaqp2 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaqp2)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *OFFSET,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *JPVT,
                    DOUBLE           *TAU,
                    DOUBLE           *VN1,
                    DOUBLE           *VN2,
                    DOUBLE           *WORK);

//-- dlaqps --------------------------------------------------------------------
void
LAPACK_IMPL(dlaqps)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *OFFSET,
                    const INTEGER    *NB,
                    INTEGER          *KB,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *JPVT,
                    DOUBLE           *TAU,
                    DOUBLE           *VN1,
                    DOUBLE           *VN2,
                    DOUBLE           *AUXV,
                    DOUBLE           *F,
                    const INTEGER    *LDF);

//-- dlaqr0 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaqr0)(const LOGICAL    *WANTT,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE           *H,
                    const INTEGER    *LDH,
                    DOUBLE           *WR,
                    DOUBLE           *WI,
                    const INTEGER    *ILOZ,
                    const INTEGER    *IHIZ,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dlaqr1 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaqr1)(const INTEGER    *N,
                    const DOUBLE     *H,
                    const INTEGER    *LDH,
                    const DOUBLE     *SR1,
                    const DOUBLE     *SI1,
                    const DOUBLE     *SR2,
                    const DOUBLE     *SI2,
                    DOUBLE           *V);

//-- dlaqr2 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaqr2)(const LOGICAL    *WANTT,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *N,
                    const INTEGER    *KTOP,
                    const INTEGER    *KBOT,
                    const INTEGER    *NW,
                    DOUBLE           *H,
                    const INTEGER    *LDH,
                    const INTEGER    *ILOZ,
                    const INTEGER    *IHIZ,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *NS,
                    INTEGER          *ND,
                    DOUBLE           *SR,
                    DOUBLE           *SI,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    const INTEGER    *NH,
                    DOUBLE           *T,
                    const INTEGER    *LDT,
                    const INTEGER    *NV,
                    DOUBLE           *WV,
                    const INTEGER    *LDWV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK);

//-- dlaqr3 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaqr3)(const LOGICAL    *WANTT,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *N,
                    const INTEGER    *KTOP,
                    const INTEGER    *KBOT,
                    const INTEGER    *NW,
                    DOUBLE           *H,
                    const INTEGER    *LDH,
                    const INTEGER    *ILOZ,
                    const INTEGER    *IHIZ,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *NS,
                    INTEGER          *ND,
                    DOUBLE           *SR,
                    DOUBLE           *SI,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    const INTEGER    *NH,
                    DOUBLE           *T,
                    const INTEGER    *LDT,
                    const INTEGER    *NV,
                    DOUBLE           *WV,
                    const INTEGER    *LDWV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK);

//-- dlaqr4 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaqr4)(const LOGICAL    *WANTT,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE           *H,
                    const INTEGER    *LDH,
                    DOUBLE           *WR,
                    DOUBLE           *WI,
                    const INTEGER    *ILOZ,
                    const INTEGER    *IHIZ,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dlaqr5 --------------------------------------------------------------------
void
LAPACK_IMPL(dlaqr5)(const LOGICAL    *WANTT,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *KACC22,
                    const INTEGER    *N,
                    const INTEGER    *KTOP,
                    const INTEGER    *KBOT,
                    const INTEGER    *NSHFTS,
                    DOUBLE           *SR,
                    DOUBLE           *SI,
                    DOUBLE           *H,
                    const INTEGER    *LDH,
                    const INTEGER    *ILOZ,
                    const INTEGER    *IHIZ,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    const INTEGER    *NV,
                    DOUBLE           *WV,
                    const INTEGER    *LDWV,
                    const INTEGER    *NH,
                    DOUBLE           *WH,
                    const INTEGER    *LDWH);

//-- dlaqsb --------------------------------------------------------------------
void
LAPACK_IMPL(dlaqsb)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- dlaqsp --------------------------------------------------------------------
void
LAPACK_IMPL(dlaqsp)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- dlaqsy --------------------------------------------------------------------
void
LAPACK_IMPL(dlaqsy)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- dlaqtr --------------------------------------------------------------------
void
LAPACK_IMPL(dlaqtr)(const LOGICAL    *LTRAN,
                    const LOGICAL    *LREAL,
                    const INTEGER    *N,
                    const DOUBLE     *T,
                    const INTEGER    *LDT,
                    const DOUBLE     *B,
                    const DOUBLE     *W,
                    DOUBLE           *SCALE,
                    DOUBLE           *X,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dlar1v --------------------------------------------------------------------
void
LAPACK_IMPL(dlar1v)(const INTEGER    *N,
                    const INTEGER    *B1,
                    const INTEGER    *BN,
                    const DOUBLE     *LAMBDA,
                    const DOUBLE     *D,
                    const DOUBLE     *L,
                    const DOUBLE     *LD,
                    const DOUBLE     *LLD,
                    const DOUBLE     *PIVMIN,
                    const DOUBLE     *GAPTOL,
                    DOUBLE           *Z,
                    const LOGICAL    *WANTNC,
                    INTEGER          *NEGCNT,
                    DOUBLE           *ZTZ,
                    DOUBLE           *MINGMA,
                    INTEGER          *R,
                    INTEGER          *ISUPPZ,
                    DOUBLE           *NRMINV,
                    DOUBLE           *RESID,
                    DOUBLE           *RQCORR,
                    DOUBLE           *WORK);

//-- dlar2v --------------------------------------------------------------------
void
LAPACK_IMPL(dlar2v)(const INTEGER    *N,
                    DOUBLE           *X,
                    DOUBLE           *Y,
                    DOUBLE           *Z,
                    const INTEGER    *INCX,
                    const DOUBLE     *C,
                    const DOUBLE     *S,
                    const INTEGER    *INCC);

//-- dlarf ---------------------------------------------------------------------
void
LAPACK_IMPL(dlarf)(const char           *SIDE,
                   const INTEGER        *M,
                   const INTEGER        *N,
                   const DOUBLE         *V,
                   const INTEGER        *INCV,
                   const DOUBLE         *TAU,
                   DOUBLE               *C,
                   const INTEGER        *LDC,
                   DOUBLE               *WORK);

//-- dlarfb --------------------------------------------------------------------
void
LAPACK_IMPL(dlarfb)(const char       *SIDE,
                    const char       *TRANS,
                    const char       *DIRECT,
                    const char       *STOREV,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const DOUBLE     *V,
                    const INTEGER    *LDV,
                    const DOUBLE     *T,
                    const INTEGER    *LDT,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    const INTEGER    *LDWORK);

//-- dlarfg --------------------------------------------------------------------
void
LAPACK_IMPL(dlarfg)(const INTEGER    *N,
                    DOUBLE           *ALPHA,
                    DOUBLE           *X,
                    const INTEGER    *INCX,
                    DOUBLE           *TAU);

//-- dlarfgp -------------------------------------------------------------------
void
LAPACK_IMPL(dlarfgp)(const INTEGER    *N,
                     DOUBLE           *ALPHA,
                     DOUBLE           *X,
                     const INTEGER    *INCX,
                     DOUBLE           *TAU);

//-- dlarft --------------------------------------------------------------------
void
LAPACK_IMPL(dlarft)(const char       *DIRECT,
                    const char       *STOREV,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    const DOUBLE     *TAU,
                    DOUBLE           *T,
                    const INTEGER    *LDT);

//-- dlarfx --------------------------------------------------------------------
void
LAPACK_IMPL(dlarfx)(const char       *SIDE,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *V,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK);

//-- dlargv --------------------------------------------------------------------
void
LAPACK_IMPL(dlargv)(const INTEGER    *N,
                    DOUBLE           *X,
                    const INTEGER    *INCX,
                    DOUBLE           *Y,
                    const INTEGER    *INCY,
                    DOUBLE           *C,
                    const INTEGER    *INCC);

//-- dlarnv --------------------------------------------------------------------
void
LAPACK_IMPL(dlarnv)(const INTEGER    *IDIST,
                    INTEGER          *ISEED,
                    const INTEGER    *N,
                    DOUBLE           *X);

//-- dlarra --------------------------------------------------------------------
void
LAPACK_IMPL(dlarra)(const INTEGER    *N,
                    const DOUBLE     *D,
                    DOUBLE           *E,
                    DOUBLE           *E2,
                    const DOUBLE     *SPLTOL,
                    const DOUBLE     *TNRM,
                    INTEGER          *NSPLIT,
                    INTEGER          *ISPLIT,
                    INTEGER          *INFO);

//-- dlarrb --------------------------------------------------------------------
void
LAPACK_IMPL(dlarrb)(const INTEGER    *N,
                    const DOUBLE     *D,
                    const DOUBLE     *LLD,
                    const INTEGER    *IFIRST,
                    const INTEGER    *ILAST,
                    const DOUBLE     *RTOL1,
                    const DOUBLE     *RTOL2,
                    const INTEGER    *OFFSET,
                    DOUBLE           *W,
                    DOUBLE           *WGAP,
                    DOUBLE           *WERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    const DOUBLE     *PIVMIN,
                    const DOUBLE     *SPDIAM,
                    const INTEGER    *TWIST,
                    INTEGER          *INFO);

//-- dlarrc --------------------------------------------------------------------
void
LAPACK_IMPL(dlarrc)(const char       *JOBT,
                    const INTEGER    *N,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    const DOUBLE     *PIVMIN,
                    INTEGER          *EIGCNT,
                    INTEGER          *LCNT,
                    INTEGER          *RCNT,
                    INTEGER          *INFO);

//-- dlarrd --------------------------------------------------------------------
void
LAPACK_IMPL(dlarrd)(const char       *RANGE,
                    const char       *ORDER,
                    const INTEGER    *N,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *GERS,
                    const DOUBLE     *RELTOL,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    const DOUBLE     *E2,
                    const DOUBLE     *PIVMIN,
                    const INTEGER    *NSPLIT,
                    const INTEGER    *ISPLIT,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE           *WERR,
                    DOUBLE           *WL,
                    DOUBLE           *WU,
                    INTEGER          *IBLOCK,
                    INTEGER          *INDEXW,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dlarre --------------------------------------------------------------------
void
LAPACK_IMPL(dlarre)(const char       *RANGE,
                    const INTEGER    *N,
                    DOUBLE           *VL,
                    DOUBLE           *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *E2,
                    const DOUBLE     *RTOL1,
                    const DOUBLE     *RTOL2,
                    const DOUBLE     *SPLTOL,
                    INTEGER          *NSPLIT,
                    INTEGER          *ISPLIT,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE           *WERR,
                    DOUBLE           *WGAP,
                    INTEGER          *IBLOCK,
                    INTEGER          *INDEXW,
                    DOUBLE           *GERS,
                    DOUBLE           *PIVMIN,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dlarrf --------------------------------------------------------------------
void
LAPACK_IMPL(dlarrf)(const INTEGER    *N,
                    const DOUBLE     *D,
                    const DOUBLE     *L,
                    const DOUBLE     *LD,
                    const INTEGER    *CLSTRT,
                    const INTEGER    *CLEND,
                    const DOUBLE     *W,
                    DOUBLE           *WGAP,
                    const DOUBLE     *WERR,
                    const DOUBLE     *SPDIAM,
                    const DOUBLE     *CLGAPL,
                    const DOUBLE     *CLGAPR,
                    const DOUBLE     *PIVMIN,
                    DOUBLE           *SIGMA,
                    DOUBLE           *DPLUS,
                    DOUBLE           *LPLUS,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dlarrj --------------------------------------------------------------------
void
LAPACK_IMPL(dlarrj)(const INTEGER    *N,
                    const DOUBLE     *D,
                    const DOUBLE     *E2,
                    const INTEGER    *IFIRST,
                    const INTEGER    *ILAST,
                    const DOUBLE     *RTOL,
                    const INTEGER    *OFFSET,
                    DOUBLE           *W,
                    DOUBLE           *WERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    const DOUBLE     *PIVMIN,
                    const DOUBLE     *SPDIAM,
                    INTEGER          *INFO);

//-- dlarrk --------------------------------------------------------------------
void
LAPACK_IMPL(dlarrk)(const INTEGER    *N,
                    const INTEGER    *IW,
                    const DOUBLE     *GL,
                    const DOUBLE     *GU,
                    const DOUBLE     *D,
                    const DOUBLE     *E2,
                    const DOUBLE     *PIVMIN,
                    const DOUBLE     *RELTOL,
                    DOUBLE           *W,
                    DOUBLE           *WERR,
                    INTEGER          *INFO);

//-- dlarrr --------------------------------------------------------------------
void
LAPACK_IMPL(dlarrr)(const INTEGER    *N,
                    const DOUBLE     *D,
                    DOUBLE           *E,
                    INTEGER          *INFO);

//-- dlarrv --------------------------------------------------------------------
void
LAPACK_IMPL(dlarrv)(const INTEGER    *N,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    DOUBLE           *D,
                    DOUBLE           *L,
                    const DOUBLE     *PIVMIN,
                    const INTEGER    *ISPLIT,
                    const INTEGER    *M,
                    const INTEGER    *DOL,
                    const INTEGER    *DOU,
                    const DOUBLE     *MINRGP,
                    const DOUBLE     *RTOL1,
                    const DOUBLE     *RTOL2,
                    DOUBLE           *W,
                    DOUBLE           *WERR,
                    DOUBLE           *WGAP,
                    const INTEGER    *IBLOCK,
                    const INTEGER    *INDEXW,
                    const DOUBLE     *GERS,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *ISUPPZ,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dlarscl2 ------------------------------------------------------------------
void
LAPACK_IMPL(dlarscl2)(const INTEGER    *M,
                      const INTEGER    *N,
                      const DOUBLE     *D,
                      DOUBLE           *X,
                      const INTEGER    *LDX);

//-- dlartg --------------------------------------------------------------------
void
LAPACK_IMPL(dlartg)(const DOUBLE     *F,
                    const DOUBLE     *G,
                    DOUBLE           *CS,
                    DOUBLE           *SN,
                    DOUBLE           *R);

//-- dlartgp -------------------------------------------------------------------
void
LAPACK_IMPL(dlartgp)(const DOUBLE     *F,
                     const DOUBLE     *G,
                     DOUBLE           *CS,
                     DOUBLE           *SN,
                     DOUBLE           *R);

//-- dlartgs -------------------------------------------------------------------
void
LAPACK_IMPL(dlartgs)(const DOUBLE     *X,
                     const DOUBLE     *Y,
                     const DOUBLE     *SIGMA,
                     DOUBLE           *CS,
                     DOUBLE           *SN);

//-- dlartv --------------------------------------------------------------------
void
LAPACK_IMPL(dlartv)(const INTEGER    *N,
                    DOUBLE           *X,
                    const INTEGER    *INCX,
                    DOUBLE           *Y,
                    const INTEGER    *INCY,
                    const DOUBLE     *C,
                    const DOUBLE     *S,
                    const INTEGER    *INCC);

//-- dlaruv --------------------------------------------------------------------
void
LAPACK_IMPL(dlaruv)(INTEGER          *ISEED,
                    const INTEGER    *N,
                    DOUBLE           *X);

//-- dlarz ---------------------------------------------------------------------
void
LAPACK_IMPL(dlarz)(const char           *SIDE,
                   const INTEGER        *M,
                   const INTEGER        *N,
                   const INTEGER        *L,
                   const DOUBLE         *V,
                   const INTEGER        *INCV,
                   const DOUBLE         *TAU,
                   DOUBLE               *C,
                   const INTEGER        *LDC,
                   DOUBLE               *WORK);

//-- dlarzb --------------------------------------------------------------------
void
LAPACK_IMPL(dlarzb)(const char       *SIDE,
                    const char       *TRANS,
                    const char       *DIRECT,
                    const char       *STOREV,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const INTEGER    *L,
                    const DOUBLE     *V,
                    const INTEGER    *LDV,
                    const DOUBLE     *T,
                    const INTEGER    *LDT,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    const INTEGER    *LDWORK);

//-- dlarzt --------------------------------------------------------------------
void
LAPACK_IMPL(dlarzt)(const char       *DIRECT,
                    const char       *STOREV,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    const DOUBLE     *TAU,
                    DOUBLE           *T,
                    const INTEGER    *LDT);

//-- dlas2 ---------------------------------------------------------------------
void
LAPACK_IMPL(dlas2)(const DOUBLE     *F,
                   const DOUBLE     *G,
                   const DOUBLE     *H,
                   DOUBLE           *SSMIN,
                   DOUBLE           *SSMAX);

//-- dlascl --------------------------------------------------------------------
void
LAPACK_IMPL(dlascl)(const char       *TYPE,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    const DOUBLE     *CFROM,
                    const DOUBLE     *CTO,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dlascl2 -------------------------------------------------------------------
void
LAPACK_IMPL(dlascl2)(const INTEGER    *M,
                     const INTEGER    *N,
                     const DOUBLE     *D,
                     DOUBLE           *X,
                     const INTEGER    *LDX);

//-- dlasd0 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasd0)(const INTEGER    *N,
                    const INTEGER    *SQRE,
                    DOUBLE           *D,
                    const DOUBLE     *E,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    DOUBLE           *VT,
                    const INTEGER    *LDVT,
                    const INTEGER    *SMLSIZ,
                    INTEGER          *IWORK,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dlasd1 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasd1)(const INTEGER    *NL,
                    const INTEGER    *NR,
                    const INTEGER    *SQRE,
                    DOUBLE           *D,
                    DOUBLE           *ALPHA,
                    DOUBLE           *BETA,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    DOUBLE           *VT,
                    const INTEGER    *LDVT,
                    INTEGER          *IDXQ,
                    INTEGER          *IWORK,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dlasd2 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasd2)(const INTEGER    *NL,
                    const INTEGER    *NR,
                    const INTEGER    *SQRE,
                    INTEGER          *K,
                    DOUBLE           *D,
                    DOUBLE           *Z,
                    const DOUBLE     *ALPHA,
                    const DOUBLE     *BETA,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    DOUBLE           *VT,
                    const INTEGER    *LDVT,
                    DOUBLE           *DSIGMA,
                    DOUBLE           *U2,
                    const INTEGER    *LDU2,
                    DOUBLE           *VT2,
                    const INTEGER    *LDVT2,
                    INTEGER          *IDXP,
                    INTEGER          *IDX,
                    INTEGER          *IDXC,
                    INTEGER          *IDXQ,
                    INTEGER          *COLTYP,
                    INTEGER          *INFO);

//-- dlasd3 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasd3)(const INTEGER    *NL,
                    const INTEGER    *NR,
                    const INTEGER    *SQRE,
                    const INTEGER    *K,
                    DOUBLE           *D,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    const DOUBLE     *DSIGMA,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    DOUBLE           *U2,
                    const INTEGER    *LDU2,
                    DOUBLE           *VT,
                    const INTEGER    *LDVT,
                    DOUBLE           *VT2,
                    const INTEGER    *LDVT2,
                    const INTEGER    *IDXC,
                    const INTEGER    *CTOT,
                    const DOUBLE     *Z,
                    INTEGER          *INFO);

//-- dlasd4 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasd4)(const INTEGER    *N,
                    const INTEGER    *I,
                    const DOUBLE     *D,
                    const DOUBLE     *Z,
                    DOUBLE           *DELTA,
                    const DOUBLE     *RHO,
                    DOUBLE           *SIGMA,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dlasd5 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasd5)(const INTEGER    *I,
                    const DOUBLE     *D,
                    const DOUBLE     *Z,
                    DOUBLE           *DELTA,
                    const DOUBLE     *RHO,
                    DOUBLE           *DSIGMA,
                    DOUBLE           *WORK);

//-- dlasd6 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasd6)(const INTEGER    *ICOMPQ,
                    const INTEGER    *NL,
                    const INTEGER    *NR,
                    const INTEGER    *SQRE,
                    DOUBLE           *D,
                    DOUBLE           *VF,
                    DOUBLE           *VL,
                    DOUBLE           *ALPHA,
                    DOUBLE           *BETA,
                    INTEGER          *IDXQ,
                    INTEGER          *PERM,
                    INTEGER          *GIVPTR,
                    INTEGER          *GIVCOL,
                    const INTEGER    *LDGCOL,
                    DOUBLE           *GIVNUM,
                    const INTEGER    *LDGNUM,
                    DOUBLE           *POLES,
                    DOUBLE           *DIFL,
                    DOUBLE           *DIFR,
                    DOUBLE           *Z,
                    INTEGER          *K,
                    DOUBLE           *C,
                    DOUBLE           *S,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dlasd7 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasd7)(const INTEGER    *ICOMPQ,
                    const INTEGER    *NL,
                    const INTEGER    *NR,
                    const INTEGER    *SQRE,
                    INTEGER          *K,
                    DOUBLE           *D,
                    DOUBLE           *Z,
                    DOUBLE           *ZW,
                    DOUBLE           *VF,
                    DOUBLE           *VFW,
                    DOUBLE           *VL,
                    DOUBLE           *VLW,
                    const DOUBLE     *ALPHA,
                    const DOUBLE     *BETA,
                    DOUBLE           *DSIGMA,
                    INTEGER          *IDX,
                    INTEGER          *IDXP,
                    const INTEGER    *IDXQ,
                    INTEGER          *PERM,
                    INTEGER          *GIVPTR,
                    INTEGER          *GIVCOL,
                    const INTEGER    *LDGCOL,
                    DOUBLE           *GIVNUM,
                    const INTEGER    *LDGNUM,
                    DOUBLE           *C,
                    DOUBLE           *S,
                    INTEGER          *INFO);

//-- dlasd8 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasd8)(const INTEGER    *ICOMPQ,
                    const INTEGER    *K,
                    DOUBLE           *D,
                    DOUBLE           *Z,
                    DOUBLE           *VF,
                    DOUBLE           *VL,
                    DOUBLE           *DIFL,
                    DOUBLE           *DIFR,
                    const INTEGER    *LDDIFR,
                    DOUBLE           *DSIGMA,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dlasda --------------------------------------------------------------------
void
LAPACK_IMPL(dlasda)(const INTEGER    *ICOMPQ,
                    const INTEGER    *SMLSIZ,
                    const INTEGER    *N,
                    const INTEGER    *SQRE,
                    DOUBLE           *D,
                    const DOUBLE     *E,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    DOUBLE           *VT,
                    INTEGER          *K,
                    DOUBLE           *DIFL,
                    DOUBLE           *DIFR,
                    DOUBLE           *Z,
                    DOUBLE           *POLES,
                    INTEGER          *GIVPTR,
                    INTEGER          *GIVCOL,
                    const INTEGER    *LDGCOL,
                    INTEGER          *PERM,
                    DOUBLE           *GIVNUM,
                    DOUBLE           *C,
                    DOUBLE           *S,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dlasdq --------------------------------------------------------------------
void
LAPACK_IMPL(dlasdq)(const char       *UPLO,
                    const INTEGER    *SQRE,
                    const INTEGER    *N,
                    const INTEGER    *NCVT,
                    const INTEGER    *NRU,
                    const INTEGER    *NCC,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *VT,
                    const INTEGER    *LDVT,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dlasdt --------------------------------------------------------------------
void
LAPACK_IMPL(dlasdt)(const INTEGER    *N,
                    INTEGER          *LVL,
                    INTEGER          *ND,
                    INTEGER          *INODE,
                    INTEGER          *NDIML,
                    INTEGER          *NDIMR,
                    const INTEGER    *MSUB);

//-- dlaset --------------------------------------------------------------------
void
LAPACK_IMPL(dlaset)(const char       *UPLO,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *ALPHA,
                    const DOUBLE     *BETA,
                    DOUBLE           *A,
                    const INTEGER    *LDA);

//-- dlasq1 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasq1)(const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dlasq2 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasq2)(const INTEGER    *N,
                    DOUBLE           *Z,
                    INTEGER          *INFO);

//-- dlasq3 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasq3)(const INTEGER    *I0,
                    INTEGER          *N0,
                    const DOUBLE     *Z,
                    INTEGER          *PP,
                    DOUBLE           *DMIN,
                    DOUBLE           *SIGMA,
                    DOUBLE           *DESIG,
                    const DOUBLE     *QMAX,
                    INTEGER          *NFAIL,
                    INTEGER          *ITER,
                    INTEGER          *NDIV,
                    const LOGICAL    *IEEE,
                    INTEGER          *TTYPE,
                    DOUBLE           *DMIN1,
                    DOUBLE           *DMIN2,
                    DOUBLE           *DN,
                    DOUBLE           *DN1,
                    DOUBLE           *DN2,
                    DOUBLE           *G,
                    DOUBLE           *TAU);

//-- dlasq4 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasq4)(const INTEGER    *I0,
                    const INTEGER    *N0,
                    const DOUBLE     *Z,
                    const INTEGER    *PP,
                    const INTEGER    *N0IN,
                    const DOUBLE     *DMIN,
                    const DOUBLE     *DMIN1,
                    const DOUBLE     *DMIN2,
                    const DOUBLE     *DN,
                    const DOUBLE     *DN1,
                    const DOUBLE     *DN2,
                    DOUBLE           *TAU,
                    INTEGER          *TTYPE,
                    DOUBLE           *G);

//-- dlasq5 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasq5)(const INTEGER    *I0,
                    const INTEGER    *N0,
                    const DOUBLE     *Z,
                    const INTEGER    *PP,
                    const DOUBLE     *TAU,
                    DOUBLE           *DMIN,
                    DOUBLE           *DMIN1,
                    DOUBLE           *DMIN2,
                    DOUBLE           *DN,
                    DOUBLE           *DNM1,
                    DOUBLE           *DNM2,
                    const LOGICAL    *IEEE);

//-- dlasq6 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasq6)(const INTEGER    *I0,
                    const INTEGER    *N0,
                    const DOUBLE     *Z,
                    const INTEGER    *PP,
                    DOUBLE           *DMIN,
                    DOUBLE           *DMIN1,
                    DOUBLE           *DMIN2,
                    DOUBLE           *DN,
                    DOUBLE           *DNM1,
                    DOUBLE           *DNM2);

//-- dlasr ---------------------------------------------------------------------
void
LAPACK_IMPL(dlasr)(const char           *SIDE,
                   const char           *PIVOT,
                   const char           *DIRECT,
                   const INTEGER        *M,
                   const INTEGER        *N,
                   const DOUBLE         *C,
                   const DOUBLE         *S,
                   DOUBLE               *A,
                   const INTEGER        *LDA);

//-- dlasrt --------------------------------------------------------------------
void
LAPACK_IMPL(dlasrt)(const char       *ID,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    INTEGER          *INFO);

//-- dlassq --------------------------------------------------------------------
void
LAPACK_IMPL(dlassq)(const INTEGER    *N,
                    const DOUBLE     *X,
                    const INTEGER    *INCX,
                    DOUBLE           *SCALE,
                    DOUBLE           *SUMSQ);

//-- dlasv2 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasv2)(const DOUBLE     *F,
                    const DOUBLE     *G,
                    const DOUBLE     *H,
                    DOUBLE           *SSMIN,
                    DOUBLE           *SSMAX,
                    DOUBLE           *SNR,
                    DOUBLE           *CSR,
                    DOUBLE           *SNL,
                    DOUBLE           *CSL);

//-- dlaswp --------------------------------------------------------------------
void
LAPACK_IMPL(dlaswp)(const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const INTEGER    *K1,
                    const INTEGER    *K2,
                    const INTEGER    *IPIV,
                    const INTEGER    *INCX);

//-- dlasy2 --------------------------------------------------------------------
void
LAPACK_IMPL(dlasy2)(const LOGICAL    *LTRANL,
                    const LOGICAL    *LTRANR,
                    const INTEGER    *ISGN,
                    const INTEGER    *N1,
                    const INTEGER    *N2,
                    const DOUBLE     *TL,
                    const INTEGER    *LDTL,
                    const DOUBLE     *TR,
                    const INTEGER    *LDTR,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *SCALE,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *XNORM,
                    INTEGER          *INFO);

//-- dlasyf --------------------------------------------------------------------
void
LAPACK_IMPL(dlasyf)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NB,
                    INTEGER          *KB,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    DOUBLE           *W,
                    const INTEGER    *LDW,
                    INTEGER          *INFO);

//-- dlat2s --------------------------------------------------------------------
void
LAPACK_IMPL(dlat2s)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    FLOAT            *SA,
                    const INTEGER    *LDSA,
                    INTEGER          *INFO);

//-- dlatbs --------------------------------------------------------------------
void
LAPACK_IMPL(dlatbs)(const char       *UPLO,
                    const char       *TRANS,
                    const char       *DIAG,
                    const char       *NORMIN,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *X,
                    DOUBLE           *SCALE,
                    DOUBLE           *CNORM,
                    INTEGER          *INFO);

//-- dlatdf --------------------------------------------------------------------
void
LAPACK_IMPL(dlatdf)(const INTEGER    *IJOB,
                    const INTEGER    *N,
                    const DOUBLE     *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *RHS,
                    DOUBLE           *RDSUM,
                    DOUBLE           *RDSCAL,
                    const INTEGER    *IPIV,
                    const INTEGER    *JPIV);

//-- dlatps --------------------------------------------------------------------
void
LAPACK_IMPL(dlatps)(const char       *UPLO,
                    const char       *TRANS,
                    const char       *DIAG,
                    const char       *NORMIN,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *X,
                    DOUBLE           *SCALE,
                    DOUBLE           *CNORM,
                    INTEGER          *INFO);

//-- dlatrd --------------------------------------------------------------------
void
LAPACK_IMPL(dlatrd)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NB,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *E,
                    DOUBLE           *TAU,
                    DOUBLE           *W,
                    const INTEGER    *LDW);

//-- dlatrs --------------------------------------------------------------------
void
LAPACK_IMPL(dlatrs)(const char       *UPLO,
                    const char       *TRANS,
                    const char       *DIAG,
                    const char       *NORMIN,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *X,
                    DOUBLE           *SCALE,
                    DOUBLE           *CNORM,
                    INTEGER          *INFO);

//-- dlatrz --------------------------------------------------------------------
void
LAPACK_IMPL(dlatrz)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *L,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK);

//-- dlatzm --------------------------------------------------------------------
void
LAPACK_IMPL(dlatzm)(const char       *SIDE,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *V,
                    const INTEGER    *INCV,
                    const DOUBLE     *TAU,
                    DOUBLE           *C1,
                    DOUBLE           *C2,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK);

//-- dlauu2 --------------------------------------------------------------------
void
LAPACK_IMPL(dlauu2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dlauum --------------------------------------------------------------------
void
LAPACK_IMPL(dlauum)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dopgtr --------------------------------------------------------------------
void
LAPACK_IMPL(dopgtr)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    const DOUBLE     *TAU,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dopmtr --------------------------------------------------------------------
void
LAPACK_IMPL(dopmtr)(const char       *SIDE,
                    const char       *UPLO,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dorbdb --------------------------------------------------------------------
void
LAPACK_IMPL(dorbdb)(const char       *TRANS,
                    const char       *SIGNS,
                    const INTEGER    *M,
                    const INTEGER    *P,
                    const INTEGER    *Q,
                    DOUBLE           *X11,
                    const INTEGER    *LDX11,
                    DOUBLE           *X12,
                    const INTEGER    *LDX12,
                    DOUBLE           *X21,
                    const INTEGER    *LDX21,
                    DOUBLE           *X22,
                    const INTEGER    *LDX22,
                    DOUBLE           *THETA,
                    DOUBLE           *PHI,
                    DOUBLE           *TAUP1,
                    DOUBLE           *TAUP2,
                    DOUBLE           *TAUQ1,
                    DOUBLE           *TAUQ2,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dorcsd --------------------------------------------------------------------
void
LAPACK_IMPL(dorcsd)(const char       *JOBU1,
                    const char       *JOBU2,
                    const char       *JOBV1T,
                    const char       *JOBV2T,
                    const char       *TRANS,
                    const char       *SIGNS,
                    const INTEGER    *M,
                    const INTEGER    *P,
                    const INTEGER    *Q,
                    const DOUBLE     *X11,
                    const INTEGER    *LDX11,
                    const DOUBLE     *X12,
                    const INTEGER    *LDX12,
                    const DOUBLE     *X21,
                    const INTEGER    *LDX21,
                    const DOUBLE     *X22,
                    const INTEGER    *LDX22,
                    DOUBLE           *THETA,
                    DOUBLE           *U1,
                    const INTEGER    *LDU1,
                    DOUBLE           *U2,
                    const INTEGER    *LDU2,
                    DOUBLE           *V1T,
                    const INTEGER    *LDV1T,
                    DOUBLE           *V2T,
                    const INTEGER    *LDV2T,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dorg2l --------------------------------------------------------------------
void
LAPACK_IMPL(dorg2l)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dorg2r --------------------------------------------------------------------
void
LAPACK_IMPL(dorg2r)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dorgbr --------------------------------------------------------------------
void
LAPACK_IMPL(dorgbr)(const char       *VECT,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dorghr --------------------------------------------------------------------
void
LAPACK_IMPL(dorghr)(const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dorgl2 --------------------------------------------------------------------
void
LAPACK_IMPL(dorgl2)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dorglq --------------------------------------------------------------------
void
LAPACK_IMPL(dorglq)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dorgql --------------------------------------------------------------------
void
LAPACK_IMPL(dorgql)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dorgqr --------------------------------------------------------------------
void
LAPACK_IMPL(dorgqr)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dorgr2 --------------------------------------------------------------------
void
LAPACK_IMPL(dorgr2)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dorgrq --------------------------------------------------------------------
void
LAPACK_IMPL(dorgrq)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dorgtr --------------------------------------------------------------------
void
LAPACK_IMPL(dorgtr)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dorm2l --------------------------------------------------------------------
void
LAPACK_IMPL(dorm2l)(const char       *SIDE,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dorm2r --------------------------------------------------------------------
void
LAPACK_IMPL(dorm2r)(const char       *SIDE,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dormbr --------------------------------------------------------------------
void
LAPACK_IMPL(dormbr)(const char       *VECT,
                    const char       *SIDE,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dormhr --------------------------------------------------------------------
void
LAPACK_IMPL(dormhr)(const char       *SIDE,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dorml2 --------------------------------------------------------------------
void
LAPACK_IMPL(dorml2)(const char       *SIDE,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dormlq --------------------------------------------------------------------
void
LAPACK_IMPL(dormlq)(const char       *SIDE,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dormql --------------------------------------------------------------------
void
LAPACK_IMPL(dormql)(const char       *SIDE,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dormqr --------------------------------------------------------------------
void
LAPACK_IMPL(dormqr)(const char       *SIDE,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dormr2 --------------------------------------------------------------------
void
LAPACK_IMPL(dormr2)(const char       *SIDE,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dormr3 --------------------------------------------------------------------
void
LAPACK_IMPL(dormr3)(const char       *SIDE,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const INTEGER    *L,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dormrq --------------------------------------------------------------------
void
LAPACK_IMPL(dormrq)(const char       *SIDE,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dormrz --------------------------------------------------------------------
void
LAPACK_IMPL(dormrz)(const char       *SIDE,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const INTEGER    *L,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dormtr --------------------------------------------------------------------
void
LAPACK_IMPL(dormtr)(const char       *SIDE,
                    const char       *UPLO,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dpbcon --------------------------------------------------------------------
void
LAPACK_IMPL(dpbcon)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dpbequ --------------------------------------------------------------------
void
LAPACK_IMPL(dpbequ)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *S,
                    DOUBLE           *SCOND,
                    DOUBLE           *AMAX,
                    INTEGER          *INFO);

//-- dpbrfs --------------------------------------------------------------------
void
LAPACK_IMPL(dpbrfs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    const DOUBLE     *AFB,
                    const INTEGER    *LDAFB,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dpbstf --------------------------------------------------------------------
void
LAPACK_IMPL(dpbstf)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO);

//-- dpbsv ---------------------------------------------------------------------
void
LAPACK_IMPL(dpbsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *KD,
                   const INTEGER        *NRHS,
                   DOUBLE               *AB,
                   const INTEGER        *LDAB,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- dpbsvx --------------------------------------------------------------------
void
LAPACK_IMPL(dpbsvx)(const char       *FACT,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    const INTEGER    *NRHS,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *AFB,
                    const INTEGER    *LDAFB,
                    char             *EQUED,
                    DOUBLE           *S,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *RCOND,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dpbtf2 --------------------------------------------------------------------
void
LAPACK_IMPL(dpbtf2)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO);

//-- dpbtrf --------------------------------------------------------------------
void
LAPACK_IMPL(dpbtrf)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO);

//-- dpbtrs --------------------------------------------------------------------
void
LAPACK_IMPL(dpbtrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dpftrf --------------------------------------------------------------------
void
LAPACK_IMPL(dpftrf)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    INTEGER          *INFO);

//-- dpftri --------------------------------------------------------------------
void
LAPACK_IMPL(dpftri)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    INTEGER          *INFO);

//-- dpftrs --------------------------------------------------------------------
void
LAPACK_IMPL(dpftrs)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dpocon --------------------------------------------------------------------
void
LAPACK_IMPL(dpocon)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dpoequ --------------------------------------------------------------------
void
LAPACK_IMPL(dpoequ)(const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *S,
                    DOUBLE           *SCOND,
                    DOUBLE           *AMAX,
                    INTEGER          *INFO);

//-- dpoequb -------------------------------------------------------------------
void
LAPACK_IMPL(dpoequb)(const INTEGER    *N,
                     const DOUBLE     *A,
                     const INTEGER    *LDA,
                     DOUBLE           *S,
                     DOUBLE           *SCOND,
                     DOUBLE           *AMAX,
                     INTEGER          *INFO);

//-- dporfs --------------------------------------------------------------------
void
LAPACK_IMPL(dporfs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *AF,
                    const INTEGER    *LDAF,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dposv ---------------------------------------------------------------------
void
LAPACK_IMPL(dposv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- dposvx --------------------------------------------------------------------
void
LAPACK_IMPL(dposvx)(const char       *FACT,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *AF,
                    const INTEGER    *LDAF,
                    char             *EQUED,
                    DOUBLE           *S,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *RCOND,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dpotf2 --------------------------------------------------------------------
void
LAPACK_IMPL(dpotf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dpotrf --------------------------------------------------------------------
void
LAPACK_IMPL(dpotrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dpotri --------------------------------------------------------------------
void
LAPACK_IMPL(dpotri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dpotrs --------------------------------------------------------------------
void
LAPACK_IMPL(dpotrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dppcon --------------------------------------------------------------------
void
LAPACK_IMPL(dppcon)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dppequ --------------------------------------------------------------------
void
LAPACK_IMPL(dppequ)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *S,
                    DOUBLE           *SCOND,
                    DOUBLE           *AMAX,
                    INTEGER          *INFO);

//-- dpprfs --------------------------------------------------------------------
void
LAPACK_IMPL(dpprfs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AP,
                    const DOUBLE     *AFP,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dppsv ---------------------------------------------------------------------
void
LAPACK_IMPL(dppsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *AP,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- dppsvx --------------------------------------------------------------------
void
LAPACK_IMPL(dppsvx)(const char       *FACT,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE           *AP,
                    DOUBLE           *AFP,
                    char             *EQUED,
                    DOUBLE           *S,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *RCOND,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dpptrf --------------------------------------------------------------------
void
LAPACK_IMPL(dpptrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    INTEGER          *INFO);

//-- dpptri --------------------------------------------------------------------
void
LAPACK_IMPL(dpptri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    INTEGER          *INFO);

//-- dpptrs --------------------------------------------------------------------
void
LAPACK_IMPL(dpptrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AP,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dpstf2 --------------------------------------------------------------------
void
LAPACK_IMPL(dpstf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *PIV,
                    INTEGER          *RANK,
                    const DOUBLE     *TOL,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dpstrf --------------------------------------------------------------------
void
LAPACK_IMPL(dpstrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *PIV,
                    INTEGER          *RANK,
                    const DOUBLE     *TOL,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dptcon --------------------------------------------------------------------
void
LAPACK_IMPL(dptcon)(const INTEGER    *N,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dpteqr --------------------------------------------------------------------
void
LAPACK_IMPL(dpteqr)(const char       *COMPZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dptrfs --------------------------------------------------------------------
void
LAPACK_IMPL(dptrfs)(const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    const DOUBLE     *DF,
                    const DOUBLE     *EF,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dptsv ---------------------------------------------------------------------
void
LAPACK_IMPL(dptsv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *D,
                   DOUBLE               *E,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- dptsvx --------------------------------------------------------------------
void
LAPACK_IMPL(dptsvx)(const char       *FACT,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    DOUBLE           *DF,
                    DOUBLE           *EF,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *RCOND,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dpttrf --------------------------------------------------------------------
void
LAPACK_IMPL(dpttrf)(const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    INTEGER          *INFO);

//-- dpttrs --------------------------------------------------------------------
void
LAPACK_IMPL(dpttrs)(const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dptts2 --------------------------------------------------------------------
void
LAPACK_IMPL(dptts2)(const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    DOUBLE           *B,
                    const INTEGER    *LDB);

//-- drscl ---------------------------------------------------------------------
void
LAPACK_IMPL(drscl)(const INTEGER        *N,
                   const DOUBLE         *SA,
                   DOUBLE               *SX,
                   const INTEGER        *INCX);

//-- dsbev ---------------------------------------------------------------------
void
LAPACK_IMPL(dsbev)(const char           *JOBZ,
                   const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *KD,
                   DOUBLE               *AB,
                   const INTEGER        *LDAB,
                   DOUBLE               *W,
                   DOUBLE               *Z,
                   const INTEGER        *LDZ,
                   DOUBLE               *WORK,
                   INTEGER              *INFO);

//-- dsbevd --------------------------------------------------------------------
void
LAPACK_IMPL(dsbevd)(const char       *JOBZ,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *W,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- dsbevx --------------------------------------------------------------------
void
LAPACK_IMPL(dsbevx)(const char       *JOBZ,
                    const char       *RANGE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *IFAIL,
                    INTEGER          *INFO);

//-- dsbgst --------------------------------------------------------------------
void
LAPACK_IMPL(dsbgst)(const char       *VECT,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KA,
                    const INTEGER    *KB,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    const DOUBLE     *BB,
                    const INTEGER    *LDBB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dsbgv ---------------------------------------------------------------------
void
LAPACK_IMPL(dsbgv)(const char           *JOBZ,
                   const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *KA,
                   const INTEGER        *KB,
                   DOUBLE               *AB,
                   const INTEGER        *LDAB,
                   DOUBLE               *BB,
                   const INTEGER        *LDBB,
                   DOUBLE               *W,
                   DOUBLE               *Z,
                   const INTEGER        *LDZ,
                   DOUBLE               *WORK,
                   INTEGER              *INFO);

//-- dsbgvd --------------------------------------------------------------------
void
LAPACK_IMPL(dsbgvd)(const char       *JOBZ,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KA,
                    const INTEGER    *KB,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *BB,
                    const INTEGER    *LDBB,
                    DOUBLE           *W,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- dsbgvx --------------------------------------------------------------------
void
LAPACK_IMPL(dsbgvx)(const char       *JOBZ,
                    const char       *RANGE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KA,
                    const INTEGER    *KB,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *BB,
                    const INTEGER    *LDBB,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *IFAIL,
                    INTEGER          *INFO);

//-- dsbtrd --------------------------------------------------------------------
void
LAPACK_IMPL(dsbtrd)(const char       *VECT,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dsecnd --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dsecnd)();

//-- dsfrk ---------------------------------------------------------------------
void
LAPACK_IMPL(dsfrk)(const char           *TRANSR,
                   const char           *UPLO,
                   const char           *TRANS,
                   const INTEGER        *N,
                   const INTEGER        *K,
                   const DOUBLE         *ALPHA,
                   const DOUBLE         *A,
                   const INTEGER        *LDA,
                   const DOUBLE         *BETA,
                   DOUBLE               *C);

//-- dsgesv --------------------------------------------------------------------
void
LAPACK_IMPL(dsgesv)(const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *WORK,
                    FLOAT            *SWORK,
                    INTEGER          *ITER,
                    INTEGER          *INFO);

//-- dspcon --------------------------------------------------------------------
void
LAPACK_IMPL(dspcon)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    const INTEGER    *IPIV,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dspev ---------------------------------------------------------------------
void
LAPACK_IMPL(dspev)(const char           *JOBZ,
                   const char           *UPLO,
                   const INTEGER        *N,
                   DOUBLE               *AP,
                   DOUBLE               *W,
                   DOUBLE               *Z,
                   const INTEGER        *LDZ,
                   DOUBLE               *WORK,
                   INTEGER              *INFO);

//-- dspevd --------------------------------------------------------------------
void
LAPACK_IMPL(dspevd)(const char       *JOBZ,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    DOUBLE           *W,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- dspevx --------------------------------------------------------------------
void
LAPACK_IMPL(dspevx)(const char       *JOBZ,
                    const char       *RANGE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *IFAIL,
                    INTEGER          *INFO);

//-- dspgst --------------------------------------------------------------------
void
LAPACK_IMPL(dspgst)(const INTEGER    *ITYPE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    const DOUBLE     *BP,
                    INTEGER          *INFO);

//-- dspgv ---------------------------------------------------------------------
void
LAPACK_IMPL(dspgv)(const INTEGER        *ITYPE,
                   const char           *JOBZ,
                   const char           *UPLO,
                   const INTEGER        *N,
                   DOUBLE               *AP,
                   DOUBLE               *BP,
                   DOUBLE               *W,
                   DOUBLE               *Z,
                   const INTEGER        *LDZ,
                   DOUBLE               *WORK,
                   INTEGER              *INFO);

//-- dspgvd --------------------------------------------------------------------
void
LAPACK_IMPL(dspgvd)(const INTEGER    *ITYPE,
                    const char       *JOBZ,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    DOUBLE           *BP,
                    DOUBLE           *W,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- dspgvx --------------------------------------------------------------------
void
LAPACK_IMPL(dspgvx)(const INTEGER    *ITYPE,
                    const char       *JOBZ,
                    const char       *RANGE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    DOUBLE           *BP,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *IFAIL,
                    INTEGER          *INFO);

//-- dsposv --------------------------------------------------------------------
void
LAPACK_IMPL(dsposv)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *WORK,
                    FLOAT            *SWORK,
                    INTEGER          *ITER,
                    INTEGER          *INFO);

//-- dsprfs --------------------------------------------------------------------
void
LAPACK_IMPL(dsprfs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AP,
                    const DOUBLE     *AFP,
                    const INTEGER    *IPIV,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dspsv ---------------------------------------------------------------------
void
LAPACK_IMPL(dspsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *AP,
                   INTEGER              *IPIV,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- dspsvx --------------------------------------------------------------------
void
LAPACK_IMPL(dspsvx)(const char       *FACT,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AP,
                    DOUBLE           *AFP,
                    INTEGER          *IPIV,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *RCOND,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dsptrd --------------------------------------------------------------------
void
LAPACK_IMPL(dsptrd)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *TAU,
                    INTEGER          *INFO);

//-- dsptrf --------------------------------------------------------------------
void
LAPACK_IMPL(dsptrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- dsptri --------------------------------------------------------------------
void
LAPACK_IMPL(dsptri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    const INTEGER    *IPIV,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dsptrs --------------------------------------------------------------------
void
LAPACK_IMPL(dsptrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AP,
                    const INTEGER    *IPIV,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dstebz --------------------------------------------------------------------
void
LAPACK_IMPL(dstebz)(const char       *RANGE,
                    const char       *ORDER,
                    const INTEGER    *N,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    INTEGER          *M,
                    INTEGER          *NSPLIT,
                    DOUBLE           *W,
                    INTEGER          *IBLOCK,
                    INTEGER          *ISPLIT,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dstedc --------------------------------------------------------------------
void
LAPACK_IMPL(dstedc)(const char       *COMPZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- dstegr --------------------------------------------------------------------
void
LAPACK_IMPL(dstegr)(const char       *JOBZ,
                    const char       *RANGE,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *ISUPPZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- dstein --------------------------------------------------------------------
void
LAPACK_IMPL(dstein)(const INTEGER    *N,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    const INTEGER    *M,
                    const DOUBLE     *W,
                    const INTEGER    *IBLOCK,
                    const INTEGER    *ISPLIT,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *IFAIL,
                    INTEGER          *INFO);

//-- dstemr --------------------------------------------------------------------
void
LAPACK_IMPL(dstemr)(const char       *JOBZ,
                    const char       *RANGE,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    const INTEGER    *NZC,
                    INTEGER          *ISUPPZ,
                    LOGICAL          *TRYRAC,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- dsteqr --------------------------------------------------------------------
void
LAPACK_IMPL(dsteqr)(const char       *COMPZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dsterf --------------------------------------------------------------------
void
LAPACK_IMPL(dsterf)(const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    INTEGER          *INFO);

//-- dstev ---------------------------------------------------------------------
void
LAPACK_IMPL(dstev)(const char           *JOBZ,
                   const INTEGER        *N,
                   DOUBLE               *D,
                   DOUBLE               *E,
                   DOUBLE               *Z,
                   const INTEGER        *LDZ,
                   DOUBLE               *WORK,
                   INTEGER              *INFO);

//-- dstevd --------------------------------------------------------------------
void
LAPACK_IMPL(dstevd)(const char       *JOBZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- dstevr --------------------------------------------------------------------
void
LAPACK_IMPL(dstevr)(const char       *JOBZ,
                    const char       *RANGE,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *ISUPPZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- dstevx --------------------------------------------------------------------
void
LAPACK_IMPL(dstevx)(const char       *JOBZ,
                    const char       *RANGE,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *IFAIL,
                    INTEGER          *INFO);

//-- dsycon --------------------------------------------------------------------
void
LAPACK_IMPL(dsycon)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dsyconv -------------------------------------------------------------------
void
LAPACK_IMPL(dsyconv)(const char       *UPLO,
                     const char       *WAY,
                     const INTEGER    *N,
                     const DOUBLE     *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE           *WORK,
                     INTEGER          *INFO);

//-- dsyequb -------------------------------------------------------------------
void
LAPACK_IMPL(dsyequb)(const char       *UPLO,
                     const INTEGER    *N,
                     const DOUBLE     *A,
                     const INTEGER    *LDA,
                     DOUBLE           *S,
                     DOUBLE           *SCOND,
                     DOUBLE           *AMAX,
                     DOUBLE           *WORK,
                     INTEGER          *INFO);

//-- dsyev ---------------------------------------------------------------------
void
LAPACK_IMPL(dsyev)(const char           *JOBZ,
                   const char           *UPLO,
                   const INTEGER        *N,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   DOUBLE               *W,
                   DOUBLE               *WORK,
                   const INTEGER        *LWORK,
                   INTEGER              *INFO);

//-- dsyevd --------------------------------------------------------------------
void
LAPACK_IMPL(dsyevd)(const char       *JOBZ,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *W,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- dsyevr --------------------------------------------------------------------
void
LAPACK_IMPL(dsyevr)(const char       *JOBZ,
                    const char       *RANGE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *ISUPPZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- dsyevx --------------------------------------------------------------------
void
LAPACK_IMPL(dsyevx)(const char       *JOBZ,
                    const char       *RANGE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    INTEGER          *IFAIL,
                    INTEGER          *INFO);

//-- dsygs2 --------------------------------------------------------------------
void
LAPACK_IMPL(dsygs2)(const INTEGER    *ITYPE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dsygst --------------------------------------------------------------------
void
LAPACK_IMPL(dsygst)(const INTEGER    *ITYPE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dsygv ---------------------------------------------------------------------
void
LAPACK_IMPL(dsygv)(const INTEGER        *ITYPE,
                   const char           *JOBZ,
                   const char           *UPLO,
                   const INTEGER        *N,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   DOUBLE               *W,
                   DOUBLE               *WORK,
                   const INTEGER        *LWORK,
                   INTEGER              *INFO);

//-- dsygvd --------------------------------------------------------------------
void
LAPACK_IMPL(dsygvd)(const INTEGER    *ITYPE,
                    const char       *JOBZ,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *W,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- dsygvx --------------------------------------------------------------------
void
LAPACK_IMPL(dsygvx)(const INTEGER    *ITYPE,
                    const char       *JOBZ,
                    const char       *RANGE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    INTEGER          *IFAIL,
                    INTEGER          *INFO);

//-- dsyrfs --------------------------------------------------------------------
void
LAPACK_IMPL(dsyrfs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *AF,
                    const INTEGER    *LDAF,
                    const INTEGER    *IPIV,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dsysv ---------------------------------------------------------------------
void
LAPACK_IMPL(dsysv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   INTEGER              *IPIV,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   DOUBLE               *WORK,
                   const INTEGER        *LWORK,
                   INTEGER              *INFO);

//-- dsysvx --------------------------------------------------------------------
void
LAPACK_IMPL(dsysvx)(const char       *FACT,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *AF,
                    const INTEGER    *LDAF,
                    INTEGER          *IPIV,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *RCOND,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dsyswapr ------------------------------------------------------------------
void
LAPACK_IMPL(dsyswapr)(const char       *UPLO,
                      const INTEGER    *N,
                      DOUBLE           *A,
                      const INTEGER    *LDA,
                      const INTEGER    *I1,
                      const INTEGER    *I2);

//-- dsytd2 --------------------------------------------------------------------
void
LAPACK_IMPL(dsytd2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *TAU,
                    INTEGER          *INFO);

//-- dsytf2 --------------------------------------------------------------------
void
LAPACK_IMPL(dsytf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- dsytrd --------------------------------------------------------------------
void
LAPACK_IMPL(dsytrd)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dsytrf --------------------------------------------------------------------
void
LAPACK_IMPL(dsytrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dsytri --------------------------------------------------------------------
void
LAPACK_IMPL(dsytri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dsytri2 -------------------------------------------------------------------
void
LAPACK_IMPL(dsytri2)(const char       *UPLO,
                     const INTEGER    *N,
                     DOUBLE           *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE           *WORK,
                     const INTEGER    *LWORK,
                     INTEGER          *INFO);

//-- dsytri2x ------------------------------------------------------------------
void
LAPACK_IMPL(dsytri2x)(const char       *UPLO,
                      const INTEGER    *N,
                      DOUBLE           *A,
                      const INTEGER    *LDA,
                      const INTEGER    *IPIV,
                      DOUBLE           *WORK,
                      const INTEGER    *NB,
                      INTEGER          *INFO);

//-- dsytrs --------------------------------------------------------------------
void
LAPACK_IMPL(dsytrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dsytrs2 -------------------------------------------------------------------
void
LAPACK_IMPL(dsytrs2)(const char       *UPLO,
                     const INTEGER    *N,
                     const INTEGER    *NRHS,
                     const DOUBLE     *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE           *B,
                     const INTEGER    *LDB,
                     DOUBLE           *WORK,
                     INTEGER          *INFO);

//-- dtbcon --------------------------------------------------------------------
void
LAPACK_IMPL(dtbcon)(const char       *NORM,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dtbrfs --------------------------------------------------------------------
void
LAPACK_IMPL(dtbrfs)(const char       *UPLO,
                    const char       *TRANS,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    const DOUBLE     *X,
                    const INTEGER    *LDX,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dtbtrs --------------------------------------------------------------------
void
LAPACK_IMPL(dtbtrs)(const char       *UPLO,
                    const char       *TRANS,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dtfsm ---------------------------------------------------------------------
void
LAPACK_IMPL(dtfsm)(const char           *TRANSR,
                   const char           *SIDE,
                   const char           *UPLO,
                   const char           *TRANS,
                   const char           *DIAG,
                   const INTEGER        *M,
                   const INTEGER        *N,
                   const DOUBLE         *ALPHA,
                   const DOUBLE         *A,
                   DOUBLE               *B,
                   const INTEGER        *LDB);

//-- dtftri --------------------------------------------------------------------
void
LAPACK_IMPL(dtftri)(const char       *TRANSR,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    INTEGER          *INFO);

//-- dtfttp --------------------------------------------------------------------
void
LAPACK_IMPL(dtfttp)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *ARF,
                    DOUBLE           *AP,
                    INTEGER          *INFO);

//-- dtfttr --------------------------------------------------------------------
void
LAPACK_IMPL(dtfttr)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *ARF,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dtgevc --------------------------------------------------------------------
void
LAPACK_IMPL(dtgevc)(const char       *SIDE,
                    const char       *HOWMNY,
                    const LOGICAL    *SELECT,
                    const INTEGER    *N,
                    const DOUBLE     *S,
                    const INTEGER    *LDS,
                    const DOUBLE     *P,
                    const INTEGER    *LDP,
                    DOUBLE           *VL,
                    const INTEGER    *LDVL,
                    DOUBLE           *VR,
                    const INTEGER    *LDVR,
                    const INTEGER    *MM,
                    INTEGER          *M,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dtgex2 --------------------------------------------------------------------
void
LAPACK_IMPL(dtgex2)(const LOGICAL    *WANTQ,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    const INTEGER    *J1,
                    const INTEGER    *N1,
                    const INTEGER    *N2,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dtgexc --------------------------------------------------------------------
void
LAPACK_IMPL(dtgexc)(const LOGICAL    *WANTQ,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *IFST,
                    INTEGER          *ILST,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dtgsen --------------------------------------------------------------------
void
LAPACK_IMPL(dtgsen)(const INTEGER    *IJOB,
                    const LOGICAL    *WANTQ,
                    const LOGICAL    *WANTZ,
                    const LOGICAL    *SELECT,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *ALPHAR,
                    DOUBLE           *ALPHAI,
                    DOUBLE           *BETA,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *M,
                    DOUBLE           *PL,
                    DOUBLE           *PR,
                    DOUBLE           *DIF,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- dtgsja --------------------------------------------------------------------
void
LAPACK_IMPL(dtgsja)(const char       *JOBU,
                    const char       *JOBV,
                    const char       *JOBQ,
                    const INTEGER    *M,
                    const INTEGER    *P,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const INTEGER    *L,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    const DOUBLE     *TOLA,
                    const DOUBLE     *TOLB,
                    DOUBLE           *ALPHA,
                    DOUBLE           *BETA,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    DOUBLE           *WORK,
                    INTEGER          *NCYCLE,
                    INTEGER          *INFO);

//-- dtgsna --------------------------------------------------------------------
void
LAPACK_IMPL(dtgsna)(const char       *JOB,
                    const char       *HOWMNY,
                    const LOGICAL    *SELECT,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    const DOUBLE     *VL,
                    const INTEGER    *LDVL,
                    const DOUBLE     *VR,
                    const INTEGER    *LDVR,
                    DOUBLE           *S,
                    DOUBLE           *DIF,
                    const INTEGER    *MM,
                    INTEGER          *M,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dtgsy2 --------------------------------------------------------------------
void
LAPACK_IMPL(dtgsy2)(const char       *TRANS,
                    const INTEGER    *IJOB,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    const DOUBLE     *D,
                    const INTEGER    *LDD,
                    const DOUBLE     *E,
                    const INTEGER    *LDE,
                    DOUBLE           *F,
                    const INTEGER    *LDF,
                    DOUBLE           *SCALE,
                    DOUBLE           *RDSUM,
                    DOUBLE           *RDSCAL,
                    INTEGER          *IWORK,
                    INTEGER          *PQ,
                    INTEGER          *INFO);

//-- dtgsyl --------------------------------------------------------------------
void
LAPACK_IMPL(dtgsyl)(const char       *TRANS,
                    const INTEGER    *IJOB,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    const DOUBLE     *D,
                    const INTEGER    *LDD,
                    const DOUBLE     *E,
                    const INTEGER    *LDE,
                    DOUBLE           *F,
                    const INTEGER    *LDF,
                    DOUBLE           *SCALE,
                    DOUBLE           *DIF,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dtpcon --------------------------------------------------------------------
void
LAPACK_IMPL(dtpcon)(const char       *NORM,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dtprfs --------------------------------------------------------------------
void
LAPACK_IMPL(dtprfs)(const char       *UPLO,
                    const char       *TRANS,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AP,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    const DOUBLE     *X,
                    const INTEGER    *LDX,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dtptri --------------------------------------------------------------------
void
LAPACK_IMPL(dtptri)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    INTEGER          *INFO);

//-- dtptrs --------------------------------------------------------------------
void
LAPACK_IMPL(dtptrs)(const char       *UPLO,
                    const char       *TRANS,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AP,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dtpttf --------------------------------------------------------------------
void
LAPACK_IMPL(dtpttf)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *ARF,
                    INTEGER          *INFO);

//-- dtpttr --------------------------------------------------------------------
void
LAPACK_IMPL(dtpttr)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dtrcon --------------------------------------------------------------------
void
LAPACK_IMPL(dtrcon)(const char       *NORM,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dtrevc --------------------------------------------------------------------
void
LAPACK_IMPL(dtrevc)(const char       *SIDE,
                    const char       *HOWMNY,
                    LOGICAL          *SELECT,
                    const INTEGER    *N,
                    const DOUBLE     *T,
                    const INTEGER    *LDT,
                    DOUBLE           *VL,
                    const INTEGER    *LDVL,
                    DOUBLE           *VR,
                    const INTEGER    *LDVR,
                    const INTEGER    *MM,
                    INTEGER          *M,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dtrexc --------------------------------------------------------------------
void
LAPACK_IMPL(dtrexc)(const char       *COMPQ,
                    const INTEGER    *N,
                    DOUBLE           *T,
                    const INTEGER    *LDT,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    INTEGER          *IFST,
                    INTEGER          *ILST,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dtrrfs --------------------------------------------------------------------
void
LAPACK_IMPL(dtrrfs)(const char       *UPLO,
                    const char       *TRANS,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    const DOUBLE     *X,
                    const INTEGER    *LDX,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dtrsen --------------------------------------------------------------------
void
LAPACK_IMPL(dtrsen)(const char       *JOB,
                    const char       *COMPQ,
                    const LOGICAL    *SELECT,
                    const INTEGER    *N,
                    DOUBLE           *T,
                    const INTEGER    *LDT,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    DOUBLE           *WR,
                    DOUBLE           *WI,
                    INTEGER          *M,
                    DOUBLE           *S,
                    DOUBLE           *SEP,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- dtrsna --------------------------------------------------------------------
void
LAPACK_IMPL(dtrsna)(const char       *JOB,
                    const char       *HOWMNY,
                    const LOGICAL    *SELECT,
                    const INTEGER    *N,
                    const DOUBLE     *T,
                    const INTEGER    *LDT,
                    const DOUBLE     *VL,
                    const INTEGER    *LDVL,
                    const DOUBLE     *VR,
                    const INTEGER    *LDVR,
                    DOUBLE           *S,
                    DOUBLE           *SEP,
                    const INTEGER    *MM,
                    INTEGER          *M,
                    DOUBLE           *WORK,
                    const INTEGER    *LDWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dtrsyl --------------------------------------------------------------------
void
LAPACK_IMPL(dtrsyl)(const char       *TRANA,
                    const char       *TRANB,
                    const INTEGER    *ISGN,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *SCALE,
                    INTEGER          *INFO);

//-- dtrti2 --------------------------------------------------------------------
void
LAPACK_IMPL(dtrti2)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dtrtri --------------------------------------------------------------------
void
LAPACK_IMPL(dtrtri)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dtrtrs --------------------------------------------------------------------
void
LAPACK_IMPL(dtrtrs)(const char       *UPLO,
                    const char       *TRANS,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dtrttf --------------------------------------------------------------------
void
LAPACK_IMPL(dtrttf)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *ARF,
                    INTEGER          *INFO);

//-- dtrttp --------------------------------------------------------------------
void
LAPACK_IMPL(dtrttp)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *AP,
                    INTEGER          *INFO);

//-- dtzrqf --------------------------------------------------------------------
void
LAPACK_IMPL(dtzrqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    INTEGER          *INFO);

//-- dtzrzf --------------------------------------------------------------------
void
LAPACK_IMPL(dtzrzf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dzsum1 --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(dzsum1)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *CX,
                    const INTEGER            *INCX);

//-- ieeeck --------------------------------------------------------------------
INTEGER
LAPACK_IMPL(ieeeck)(const INTEGER    *ISPEC,
                    const FLOAT      *ZERO,
                    const FLOAT      *ONE);

//-- iladlc --------------------------------------------------------------------
INTEGER
LAPACK_IMPL(iladlc)(const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA);

//-- iladlr --------------------------------------------------------------------
INTEGER
LAPACK_IMPL(iladlr)(const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA);

//-- ilaslc --------------------------------------------------------------------
INTEGER
LAPACK_IMPL(ilaslc)(const INTEGER    *M,
                    const INTEGER    *N,
                    const FLOAT      *A,
                    const INTEGER    *LDA);

//-- ilaslr --------------------------------------------------------------------
INTEGER
LAPACK_IMPL(ilaslr)(const INTEGER    *M,
                    const INTEGER    *N,
                    const FLOAT      *A,
                    const INTEGER    *LDA);

//-- ilatrans ------------------------------------------------------------------
INTEGER
LAPACK_IMPL(ilatrans)(const char       *TRANS);

//-- ilauplo -------------------------------------------------------------------
INTEGER
LAPACK_IMPL(ilauplo)(const char   *UPLO);

//-- ilaver --------------------------------------------------------------------
void
LAPACK_IMPL(ilaver)(INTEGER  *VERS_MAJOR,
                    INTEGER  *VERS_MINOR,
                    INTEGER  *VERS_PATCH);

//-- lsame ---------------------------------------------------------------------
LOGICAL
LAPACK_IMPL(lsame)(const char       *CA,
                   const char       *CB);

//-- lsamen --------------------------------------------------------------------
LOGICAL
LAPACK_IMPL(lsamen)(const INTEGER    *N,
                    const char       *CA,
                    const char       *CB);

//-- sgetf2 --------------------------------------------------------------------
void
LAPACK_IMPL(sgetf2)(const INTEGER    *M,
                    const INTEGER    *N,
                    FLOAT            *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- sgetrf --------------------------------------------------------------------
void
LAPACK_IMPL(sgetrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    FLOAT            *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- sgetrs --------------------------------------------------------------------
void
LAPACK_IMPL(sgetrs)(const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const FLOAT      *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    FLOAT            *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- sisnan --------------------------------------------------------------------
LOGICAL
LAPACK_IMPL(sisnan)(const FLOAT  *SIN);

//-- slag2d --------------------------------------------------------------------
void
LAPACK_IMPL(slag2d)(const INTEGER    *M,
                    const INTEGER    *N,
                    const FLOAT      *SA,
                    const INTEGER    *LDSA,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- slaisnan ------------------------------------------------------------------
LOGICAL
LAPACK_IMPL(slaisnan)(const FLOAT      *SIN1,
                      const FLOAT      *SIN2);

//-- slamch --------------------------------------------------------------------
FLOAT
LAPACK_IMPL(slamch)(const char   *CMACH);

//-- slaswp --------------------------------------------------------------------
void
LAPACK_IMPL(slaswp)(const INTEGER    *N,
                    FLOAT            *A,
                    const INTEGER    *LDA,
                    const INTEGER    *K1,
                    const INTEGER    *K2,
                    const INTEGER    *IPIV,
                    const INTEGER    *INCX);

//-- spotf2 --------------------------------------------------------------------
void
LAPACK_IMPL(spotf2)(const char       *UPLO,
                    const INTEGER    *N,
                    FLOAT            *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- spotrf --------------------------------------------------------------------
void
LAPACK_IMPL(spotrf)(const char       *UPLO,
                    const INTEGER    *N,
                    FLOAT            *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- spotrs --------------------------------------------------------------------
void
LAPACK_IMPL(spotrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const FLOAT      *A,
                    const INTEGER    *LDA,
                    FLOAT            *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- zbbcsd --------------------------------------------------------------------
void
LAPACK_IMPL(zbbcsd)(const char       *JOBU1,
                    const char       *JOBU2,
                    const char       *JOBV1T,
                    const char       *JOBV2T,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *P,
                    const INTEGER    *Q,
                    DOUBLE           *THETA,
                    DOUBLE           *PHI,
                    DOUBLE_COMPLEX   *U1,
                    const INTEGER    *LDU1,
                    DOUBLE_COMPLEX   *U2,
                    const INTEGER    *LDU2,
                    DOUBLE_COMPLEX   *V1T,
                    const INTEGER    *LDV1T,
                    DOUBLE_COMPLEX   *V2T,
                    const INTEGER    *LDV2T,
                    DOUBLE           *B11D,
                    DOUBLE           *B11E,
                    DOUBLE           *B12D,
                    DOUBLE           *B12E,
                    const DOUBLE     *B21D,
                    const DOUBLE     *B21E,
                    const DOUBLE     *B22D,
                    const DOUBLE     *B22E,
                    DOUBLE           *RWORK,
                    const INTEGER    *LRWORK,
                    INTEGER          *INFO);

//-- zbdsqr --------------------------------------------------------------------
void
LAPACK_IMPL(zbdsqr)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NCVT,
                    const INTEGER    *NRU,
                    const INTEGER    *NCC,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *VT,
                    const INTEGER    *LDVT,
                    DOUBLE_COMPLEX   *U,
                    const INTEGER    *LDU,
                    DOUBLE_COMPLEX   *C,
                    const INTEGER    *LDC,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zcgesv --------------------------------------------------------------------
void
LAPACK_IMPL(zcgesv)(const INTEGER            *N,
                    const INTEGER            *NRHS,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    INTEGER                  *IPIV,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE_COMPLEX           *WORK,
                    FLOAT_COMPLEX            *SWORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *ITER,
                    INTEGER                  *INFO);

//-- zcposv --------------------------------------------------------------------
void
LAPACK_IMPL(zcposv)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE_COMPLEX           *WORK,
                    FLOAT_COMPLEX            *SWORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *ITER,
                    INTEGER                  *INFO);

//-- zdrscl --------------------------------------------------------------------
void
LAPACK_IMPL(zdrscl)(const INTEGER    *N,
                    const DOUBLE     *SA,
                    DOUBLE_COMPLEX   *SX,
                    const INTEGER    *INCX);

//-- zgbbrd --------------------------------------------------------------------
void
LAPACK_IMPL(zgbbrd)(const char       *VECT,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *NCC,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    DOUBLE_COMPLEX   *PT,
                    const INTEGER    *LDPT,
                    DOUBLE_COMPLEX   *C,
                    const INTEGER    *LDC,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zgbcon --------------------------------------------------------------------
void
LAPACK_IMPL(zgbcon)(const char               *NORM,
                    const INTEGER            *N,
                    const INTEGER            *KL,
                    const INTEGER            *KU,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    const INTEGER            *IPIV,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zgbequ --------------------------------------------------------------------
void
LAPACK_IMPL(zgbequ)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *KL,
                    const INTEGER            *KU,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE                   *R,
                    DOUBLE                   *C,
                    DOUBLE                   *ROWCND,
                    DOUBLE                   *COLCND,
                    DOUBLE                   *AMAX,
                    INTEGER                  *INFO);

//-- zgbequb -------------------------------------------------------------------
void
LAPACK_IMPL(zgbequb)(const INTEGER            *M,
                     const INTEGER            *N,
                     const INTEGER            *KL,
                     const INTEGER            *KU,
                     const DOUBLE_COMPLEX     *AB,
                     const INTEGER            *LDAB,
                     DOUBLE                   *R,
                     DOUBLE                   *C,
                     DOUBLE                   *ROWCND,
                     DOUBLE                   *COLCND,
                     DOUBLE                   *AMAX,
                     INTEGER                  *INFO);

//-- zgbrfs --------------------------------------------------------------------
void
LAPACK_IMPL(zgbrfs)(const char               *TRANS,
                    const INTEGER            *N,
                    const INTEGER            *KL,
                    const INTEGER            *KU,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    const DOUBLE_COMPLEX     *AFB,
                    const INTEGER            *LDAFB,
                    const INTEGER            *IPIV,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zgbrfsx -------------------------------------------------------------------
/*
void
LAPACK_IMPL(zgbrfsx)(const char               *TRANS,
                     const char               *EQUED,
                     const INTEGER            *N,
                     const INTEGER            *KL,
                     const INTEGER            *KU,
                     const INTEGER            *NRHS,
                     const DOUBLE_COMPLEX     *AB,
                     const INTEGER            *LDAB,
                     const DOUBLE_COMPLEX     *AFB,
                     const INTEGER            *LDAFB,
                     const INTEGER            *IPIV,
                     DOUBLE                   *R,
                     DOUBLE                   *C,
                     const DOUBLE_COMPLEX     *B,
                     const INTEGER            *LDB,
                     DOUBLE_COMPLEX           *X,
                     const INTEGER            *LDX,
                     DOUBLE                   *RCOND,
                     DOUBLE                   *BERR,
                     const INTEGER            *N_ERR_BNDS,
                     DOUBLE                   *ERR_BNDS_NORM,
                     DOUBLE                   *ERR_BNDS_COMP,
                     const INTEGER            *NPARAMS,
                     DOUBLE                   *PARAMS,
                     DOUBLE_COMPLEX           *WORK,
                     DOUBLE                   *RWORK,
                     INTEGER                  *INFO);
*/

//-- zgbsv ---------------------------------------------------------------------
void
LAPACK_IMPL(zgbsv)(const INTEGER        *N,
                   const INTEGER        *KL,
                   const INTEGER        *KU,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *AB,
                   const INTEGER        *LDAB,
                   INTEGER              *IPIV,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zgbsvx --------------------------------------------------------------------
void
LAPACK_IMPL(zgbsvx)(const char       *FACT,
                    const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    const INTEGER    *NRHS,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    DOUBLE_COMPLEX   *AFB,
                    const INTEGER    *LDAFB,
                    INTEGER          *IPIV,
                    char             *EQUED,
                    DOUBLE           *R,
                    DOUBLE           *C,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *LDX,
                    DOUBLE           *RCOND,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zgbsvxx -------------------------------------------------------------------
/*
void
LAPACK_IMPL(zgbsvxx)(const char       *FACT,
                     const char       *TRANS,
                     const INTEGER    *N,
                     const INTEGER    *KL,
                     const INTEGER    *KU,
                     const INTEGER    *NRHS,
                     DOUBLE_COMPLEX   *AB,
                     const INTEGER    *LDAB,
                     DOUBLE_COMPLEX   *AFB,
                     const INTEGER    *LDAFB,
                     INTEGER          *IPIV,
                     char             *EQUED,
                     DOUBLE           *R,
                     DOUBLE           *C,
                     DOUBLE_COMPLEX   *B,
                     const INTEGER    *LDB,
                     DOUBLE_COMPLEX   *X,
                     const INTEGER    *LDX,
                     DOUBLE           *RCOND,
                     DOUBLE           *RPVGRW,
                     DOUBLE           *BERR,
                     const INTEGER    *N_ERR_BNDS,
                     DOUBLE           *ERR_BNDS_NORM,
                     DOUBLE           *ERR_BNDS_COMP,
                     const INTEGER    *NPARAMS,
                     DOUBLE           *PARAMS,
                     DOUBLE_COMPLEX   *WORK,
                     DOUBLE           *RWORK,
                     INTEGER          *INFO);
*/
//-- zgbtf2 --------------------------------------------------------------------
void
LAPACK_IMPL(zgbtf2)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zgbtrf --------------------------------------------------------------------
void
LAPACK_IMPL(zgbtrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zgbtrs --------------------------------------------------------------------
void
LAPACK_IMPL(zgbtrs)(const char               *TRANS,
                    const INTEGER            *N,
                    const INTEGER            *KL,
                    const INTEGER            *KU,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    const INTEGER            *IPIV,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zgebak --------------------------------------------------------------------
void
LAPACK_IMPL(zgebak)(const char       *JOB,
                    const char       *SIDE,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    const DOUBLE     *SCALE,
                    const INTEGER    *M,
                    DOUBLE_COMPLEX   *V,
                    const INTEGER    *LDV,
                    INTEGER          *INFO);

//-- zgebal --------------------------------------------------------------------
void
LAPACK_IMPL(zgebal)(const char       *JOB,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *ILO,
                    INTEGER          *IHI,
                    DOUBLE           *SCALE,
                    INTEGER          *INFO);

//-- zgebd2 --------------------------------------------------------------------
void
LAPACK_IMPL(zgebd2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *TAUQ,
                    DOUBLE_COMPLEX   *TAUP,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zgebrd --------------------------------------------------------------------
void
LAPACK_IMPL(zgebrd)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *TAUQ,
                    DOUBLE_COMPLEX   *TAUP,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zgecon --------------------------------------------------------------------
void
LAPACK_IMPL(zgecon)(const char               *NORM,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zgeequ --------------------------------------------------------------------
void
LAPACK_IMPL(zgeequ)(const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *R,
                    DOUBLE                   *C,
                    DOUBLE                   *ROWCND,
                    DOUBLE                   *COLCND,
                    DOUBLE                   *AMAX,
                    INTEGER                  *INFO);

//-- zgeequb -------------------------------------------------------------------
void
LAPACK_IMPL(zgeequb)(const INTEGER            *M,
                     const INTEGER            *N,
                     const DOUBLE_COMPLEX     *A,
                     const INTEGER            *LDA,
                     DOUBLE                   *R,
                     DOUBLE                   *C,
                     DOUBLE                   *ROWCND,
                     DOUBLE                   *COLCND,
                     DOUBLE                   *AMAX,
                     INTEGER                  *INFO);

//-- zgees ---------------------------------------------------------------------
void
LAPACK_IMPL(zgees)(const char           *JOBVS,
                   const char           *SORT,
                   LOGICAL              (*SELECT)(const DOUBLE_COMPLEX *),
                   const INTEGER        *N,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   INTEGER              *SDIM,
                   DOUBLE_COMPLEX       *W,
                   DOUBLE_COMPLEX       *VS,
                   const INTEGER        *LDVS,
                   DOUBLE_COMPLEX       *WORK,
                   const INTEGER        *LWORK,
                   DOUBLE               *RWORK,
                   LOGICAL              *BWORK,
                   INTEGER              *INFO);

//-- zgeesx --------------------------------------------------------------------
void
LAPACK_IMPL(zgeesx)(const char       *JOBVS,
                    const char       *SORT,
                    LOGICAL          (*SELECT)(const DOUBLE_COMPLEX *),
                    const char       *SENSE,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *SDIM,
                    DOUBLE_COMPLEX   *W,
                    DOUBLE_COMPLEX   *VS,
                    const INTEGER    *LDVS,
                    DOUBLE           *RCONDE,
                    DOUBLE           *RCONDV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    LOGICAL          *BWORK,
                    INTEGER          *INFO);

//-- zgeev ---------------------------------------------------------------------
void
LAPACK_IMPL(zgeev)(const char           *JOBVL,
                   const char           *JOBVR,
                   const INTEGER        *N,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   DOUBLE_COMPLEX       *W,
                   DOUBLE_COMPLEX       *VL,
                   const INTEGER        *LDVL,
                   DOUBLE_COMPLEX       *VR,
                   const INTEGER        *LDVR,
                   DOUBLE_COMPLEX       *WORK,
                   const INTEGER        *LWORK,
                   DOUBLE               *RWORK,
                   INTEGER              *INFO);

//-- zgeevx --------------------------------------------------------------------
void
LAPACK_IMPL(zgeevx)(const char       *BALANC,
                    const char       *JOBVL,
                    const char       *JOBVR,
                    const char       *SENSE,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *W,
                    DOUBLE_COMPLEX   *VL,
                    const INTEGER    *LDVL,
                    DOUBLE_COMPLEX   *VR,
                    const INTEGER    *LDVR,
                    INTEGER          *ILO,
                    INTEGER          *IHI,
                    DOUBLE           *SCALE,
                    DOUBLE           *ABNRM,
                    DOUBLE           *RCONDE,
                    DOUBLE           *RCONDV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zgegs ---------------------------------------------------------------------
void
LAPACK_IMPL(zgegs)(const char           *JOBVSL,
                   const char           *JOBVSR,
                   const INTEGER        *N,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   DOUBLE_COMPLEX       *ALPHA,
                   DOUBLE_COMPLEX       *BETA,
                   DOUBLE_COMPLEX       *VSL,
                   const INTEGER        *LDVSL,
                   DOUBLE_COMPLEX       *VSR,
                   const INTEGER        *LDVSR,
                   DOUBLE_COMPLEX       *WORK,
                   const INTEGER        *LWORK,
                   DOUBLE               *RWORK,
                   INTEGER              *INFO);

//-- zgegv ---------------------------------------------------------------------
void
LAPACK_IMPL(zgegv)(const char           *JOBVL,
                   const char           *JOBVR,
                   const INTEGER        *N,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   DOUBLE_COMPLEX       *ALPHA,
                   DOUBLE_COMPLEX       *BETA,
                   DOUBLE_COMPLEX       *VL,
                   const INTEGER        *LDVL,
                   DOUBLE_COMPLEX       *VR,
                   const INTEGER        *LDVR,
                   DOUBLE_COMPLEX       *WORK,
                   const INTEGER        *LWORK,
                   DOUBLE               *RWORK,
                   INTEGER              *INFO);

//-- zgehd2 --------------------------------------------------------------------
void
LAPACK_IMPL(zgehd2)(const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zgehrd --------------------------------------------------------------------
void
LAPACK_IMPL(zgehrd)(const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zgelq2 --------------------------------------------------------------------
void
LAPACK_IMPL(zgelq2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zgelqf --------------------------------------------------------------------
void
LAPACK_IMPL(zgelqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zgels ---------------------------------------------------------------------
void
LAPACK_IMPL(zgels)(const char           *TRANS,
                   const INTEGER        *M,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   DOUBLE_COMPLEX       *WORK,
                   const INTEGER        *LWORK,
                   INTEGER              *INFO);

//-- zgelsd --------------------------------------------------------------------
void
LAPACK_IMPL(zgelsd)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    DOUBLE                   *S,
                    const DOUBLE             *RCOND,
                    INTEGER                  *RANK,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *IWORK,
                    INTEGER                  *INFO);

//-- zgelss --------------------------------------------------------------------
void
LAPACK_IMPL(zgelss)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE           *S,
                    const DOUBLE     *RCOND,
                    INTEGER          *RANK,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zgelsx --------------------------------------------------------------------
void
LAPACK_IMPL(zgelsx)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    INTEGER          *JPVT,
                    const DOUBLE     *RCOND,
                    INTEGER          *RANK,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zgelsy --------------------------------------------------------------------
void
LAPACK_IMPL(zgelsy)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    INTEGER          *JPVT,
                    const DOUBLE     *RCOND,
                    INTEGER          *RANK,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zgeql2 --------------------------------------------------------------------
void
LAPACK_IMPL(zgeql2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zgeqlf --------------------------------------------------------------------
void
LAPACK_IMPL(zgeqlf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zgeqp3 --------------------------------------------------------------------
void
LAPACK_IMPL(zgeqp3)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *JPVT,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zgeqpf --------------------------------------------------------------------
void
LAPACK_IMPL(zgeqpf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *JPVT,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zgeqr2 --------------------------------------------------------------------
void
LAPACK_IMPL(zgeqr2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zgeqr2p -------------------------------------------------------------------
void
LAPACK_IMPL(zgeqr2p)(const INTEGER    *M,
                     const INTEGER    *N,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     DOUBLE_COMPLEX   *TAU,
                     DOUBLE_COMPLEX   *WORK,
                     INTEGER          *INFO);

//-- zgeqrf --------------------------------------------------------------------
void
LAPACK_IMPL(zgeqrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zgeqrfp -------------------------------------------------------------------
void
LAPACK_IMPL(zgeqrfp)(const INTEGER    *M,
                     const INTEGER    *N,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     DOUBLE_COMPLEX   *TAU,
                     DOUBLE_COMPLEX   *WORK,
                     const INTEGER    *LWORK,
                     INTEGER          *INFO);

//-- zgerfs --------------------------------------------------------------------
void
LAPACK_IMPL(zgerfs)(const char               *TRANS,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *AF,
                    const INTEGER            *LDAF,
                    const INTEGER            *IPIV,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zgerfsx -------------------------------------------------------------------
/*
void
LAPACK_IMPL(zgerfsx)(const char               *TRANS,
                     const char               *EQUED,
                     const INTEGER            *N,
                     const INTEGER            *NRHS,
                     const DOUBLE_COMPLEX     *A,
                     const INTEGER            *LDA,
                     const DOUBLE_COMPLEX     *AF,
                     const INTEGER            *LDAF,
                     const INTEGER            *IPIV,
                     const DOUBLE             *R,
                     const DOUBLE             *C,
                     const DOUBLE_COMPLEX     *B,
                     const INTEGER            *LDB,
                     DOUBLE_COMPLEX           *X,
                     const INTEGER            *LDX,
                     DOUBLE                   *RCOND,
                     DOUBLE                   *BERR,
                     const INTEGER            *N_ERR_BNDS,
                     DOUBLE                   *ERR_BNDS_NORM,
                     DOUBLE                   *ERR_BNDS_COMP,
                     const INTEGER            *NPARAMS,
                     DOUBLE                   *PARAMS,
                     DOUBLE_COMPLEX           *WORK,
                     DOUBLE                   *RWORK,
                     INTEGER                  *INFO);
*/

//-- zgerq2 --------------------------------------------------------------------
void
LAPACK_IMPL(zgerq2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zgerqf --------------------------------------------------------------------
void
LAPACK_IMPL(zgerqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zgesc2 --------------------------------------------------------------------
void
LAPACK_IMPL(zgesc2)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *RHS,
                    const INTEGER            *IPIV,
                    const INTEGER            *JPIV,
                    DOUBLE                   *SCALE);

//-- zgesdd --------------------------------------------------------------------
void
LAPACK_IMPL(zgesdd)(const char       *JOBZ,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE           *S,
                    DOUBLE_COMPLEX   *U,
                    const INTEGER    *LDU,
                    DOUBLE_COMPLEX   *VT,
                    const INTEGER    *LDVT,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- zgesv ---------------------------------------------------------------------
void
LAPACK_IMPL(zgesv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   INTEGER              *IPIV,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zgesvd --------------------------------------------------------------------
void
LAPACK_IMPL(zgesvd)(const char       *JOBU,
                    const char       *JOBVT,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE           *S,
                    DOUBLE_COMPLEX   *U,
                    const INTEGER    *LDU,
                    DOUBLE_COMPLEX   *VT,
                    const INTEGER    *LDVT,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zgesvx --------------------------------------------------------------------
void
LAPACK_IMPL(zgesvx)(const char       *FACT,
                    const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *AF,
                    const INTEGER    *LDAF,
                    INTEGER          *IPIV,
                    char             *EQUED,
                    DOUBLE           *R,
                    DOUBLE           *C,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *LDX,
                    DOUBLE           *RCOND,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zgesvxx -------------------------------------------------------------------
/*
void
LAPACK_IMPL(zgesvxx)(const char       *FACT,
                     const char       *TRANS,
                     const INTEGER    *N,
                     const INTEGER    *NRHS,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     DOUBLE_COMPLEX   *AF,
                     const INTEGER    *LDAF,
                     INTEGER          *IPIV,
                     char             *EQUED,
                     DOUBLE           *R,
                     DOUBLE           *C,
                     DOUBLE_COMPLEX   *B,
                     const INTEGER    *LDB,
                     DOUBLE_COMPLEX   *X,
                     const INTEGER    *LDX,
                     DOUBLE           *RCOND,
                     DOUBLE           *RPVGRW,
                     DOUBLE           *BERR,
                     const INTEGER    *N_ERR_BNDS,
                     DOUBLE           *ERR_BNDS_NORM,
                     DOUBLE           *ERR_BNDS_COMP,
                     const INTEGER    *NPARAMS,
                     DOUBLE           *PARAMS,
                     DOUBLE_COMPLEX   *WORK,
                     DOUBLE           *RWORK,
                     INTEGER          *INFO);
*/
//-- zgetc2 --------------------------------------------------------------------
void
LAPACK_IMPL(zgetc2)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *JPIV,
                    INTEGER          *INFO);

//-- zgetf2 --------------------------------------------------------------------
void
LAPACK_IMPL(zgetf2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zgetrf --------------------------------------------------------------------
void
LAPACK_IMPL(zgetrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zgetri --------------------------------------------------------------------
void
LAPACK_IMPL(zgetri)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zgetrs --------------------------------------------------------------------
void
LAPACK_IMPL(zgetrs)(const char               *TRANS,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const INTEGER            *IPIV,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zggbak --------------------------------------------------------------------
void
LAPACK_IMPL(zggbak)(const char       *JOB,
                    const char       *SIDE,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    const DOUBLE     *LSCALE,
                    const DOUBLE     *RSCALE,
                    const INTEGER    *M,
                    DOUBLE_COMPLEX   *V,
                    const INTEGER    *LDV,
                    INTEGER          *INFO);

//-- zggbal --------------------------------------------------------------------
void
LAPACK_IMPL(zggbal)(const char       *JOB,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    INTEGER          *ILO,
                    INTEGER          *IHI,
                    DOUBLE           *LSCALE,
                    DOUBLE           *RSCALE,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- zgges ---------------------------------------------------------------------
void
LAPACK_IMPL(zgges)(const char           *JOBVSL,
                   const char           *JOBVSR,
                   const char           *SORT,
                   const LOGICAL        *SELCTG,
                   const INTEGER        *N,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *SDIM,
                   DOUBLE_COMPLEX       *ALPHA,
                   DOUBLE_COMPLEX       *BETA,
                   DOUBLE_COMPLEX       *VSL,
                   const INTEGER        *LDVSL,
                   DOUBLE_COMPLEX       *VSR,
                   const INTEGER        *LDVSR,
                   DOUBLE_COMPLEX       *WORK,
                   const INTEGER        *LWORK,
                   DOUBLE               *RWORK,
                   LOGICAL              *BWORK,
                   INTEGER              *INFO);

//-- zggesx --------------------------------------------------------------------
void
LAPACK_IMPL(zggesx)(const char       *JOBVSL,
                    const char       *JOBVSR,
                    const char       *SORT,
                    const LOGICAL    *SELCTG,
                    const char       *SENSE,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    INTEGER          *SDIM,
                    DOUBLE_COMPLEX   *ALPHA,
                    DOUBLE_COMPLEX   *BETA,
                    DOUBLE_COMPLEX   *VSL,
                    const INTEGER    *LDVSL,
                    DOUBLE_COMPLEX   *VSR,
                    const INTEGER    *LDVSR,
                    DOUBLE           *RCONDE,
                    DOUBLE           *RCONDV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    LOGICAL          *BWORK,
                    INTEGER          *INFO);

//-- zggev ---------------------------------------------------------------------
void
LAPACK_IMPL(zggev)(const char           *JOBVL,
                   const char           *JOBVR,
                   const INTEGER        *N,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   DOUBLE_COMPLEX       *ALPHA,
                   DOUBLE_COMPLEX       *BETA,
                   DOUBLE_COMPLEX       *VL,
                   const INTEGER        *LDVL,
                   DOUBLE_COMPLEX       *VR,
                   const INTEGER        *LDVR,
                   DOUBLE_COMPLEX       *WORK,
                   const INTEGER        *LWORK,
                   DOUBLE               *RWORK,
                   INTEGER              *INFO);

//-- zggevx --------------------------------------------------------------------
void
LAPACK_IMPL(zggevx)(const char       *BALANC,
                    const char       *JOBVL,
                    const char       *JOBVR,
                    const char       *SENSE,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *ALPHA,
                    DOUBLE_COMPLEX   *BETA,
                    DOUBLE_COMPLEX   *VL,
                    const INTEGER    *LDVL,
                    DOUBLE_COMPLEX   *VR,
                    const INTEGER    *LDVR,
                    INTEGER          *ILO,
                    INTEGER          *IHI,
                    DOUBLE           *LSCALE,
                    DOUBLE           *RSCALE,
                    DOUBLE           *ABNRM,
                    DOUBLE           *BBNRM,
                    DOUBLE           *RCONDE,
                    DOUBLE           *RCONDV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    INTEGER          *IWORK,
                    LOGICAL          *BWORK,
                    INTEGER          *INFO);

//-- zggglm --------------------------------------------------------------------
void
LAPACK_IMPL(zggglm)(const INTEGER    *N,
                    const INTEGER    *M,
                    const INTEGER    *P,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *D,
                    DOUBLE_COMPLEX   *X,
                    DOUBLE_COMPLEX   *Y,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zgghrd --------------------------------------------------------------------
void
LAPACK_IMPL(zgghrd)(const char       *COMPQ,
                    const char       *COMPZ,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *INFO);

//-- zgglse --------------------------------------------------------------------
void
LAPACK_IMPL(zgglse)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *P,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *C,
                    DOUBLE_COMPLEX   *D,
                    DOUBLE_COMPLEX   *X,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zggqrf --------------------------------------------------------------------
void
LAPACK_IMPL(zggqrf)(const INTEGER    *N,
                    const INTEGER    *M,
                    const INTEGER    *P,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAUA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *TAUB,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zggrqf --------------------------------------------------------------------
void
LAPACK_IMPL(zggrqf)(const INTEGER    *M,
                    const INTEGER    *P,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAUA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *TAUB,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zggsvd --------------------------------------------------------------------
void
LAPACK_IMPL(zggsvd)(const char       *JOBU,
                    const char       *JOBV,
                    const char       *JOBQ,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *P,
                    INTEGER          *K,
                    INTEGER          *L,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE           *ALPHA,
                    DOUBLE           *BETA,
                    DOUBLE_COMPLEX   *U,
                    const INTEGER    *LDU,
                    DOUBLE_COMPLEX   *V,
                    const INTEGER    *LDV,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- zggsvp --------------------------------------------------------------------
void
LAPACK_IMPL(zggsvp)(const char       *JOBU,
                    const char       *JOBV,
                    const char       *JOBQ,
                    const INTEGER    *M,
                    const INTEGER    *P,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    const DOUBLE     *TOLA,
                    const DOUBLE     *TOLB,
                    INTEGER          *K,
                    INTEGER          *L,
                    DOUBLE_COMPLEX   *U,
                    const INTEGER    *LDU,
                    DOUBLE_COMPLEX   *V,
                    const INTEGER    *LDV,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    INTEGER          *IWORK,
                    DOUBLE           *RWORK,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zgtcon --------------------------------------------------------------------
void
LAPACK_IMPL(zgtcon)(const char               *NORM,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *DL,
                    const DOUBLE_COMPLEX     *D,
                    const DOUBLE_COMPLEX     *DU,
                    const DOUBLE_COMPLEX     *DU2,
                    const INTEGER            *IPIV,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zgtrfs --------------------------------------------------------------------
void
LAPACK_IMPL(zgtrfs)(const char               *TRANS,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *DL,
                    const DOUBLE_COMPLEX     *D,
                    const DOUBLE_COMPLEX     *DU,
                    const DOUBLE_COMPLEX     *DLF,
                    const DOUBLE_COMPLEX     *DF,
                    const DOUBLE_COMPLEX     *DUF,
                    const DOUBLE_COMPLEX     *DU2,
                    const INTEGER            *IPIV,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zgtsv ---------------------------------------------------------------------
void
LAPACK_IMPL(zgtsv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *DL,
                   DOUBLE_COMPLEX       *D,
                   DOUBLE_COMPLEX       *DU,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zgtsvx --------------------------------------------------------------------
void
LAPACK_IMPL(zgtsvx)(const char               *FACT,
                    const char               *TRANS,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *DL,
                    const DOUBLE_COMPLEX     *D,
                    const DOUBLE_COMPLEX     *DU,
                    DOUBLE_COMPLEX           *DLF,
                    DOUBLE_COMPLEX           *DF,
                    DOUBLE_COMPLEX           *DUF,
                    DOUBLE_COMPLEX           *DU2,
                    INTEGER                  *IPIV,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *RCOND,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zgttrf --------------------------------------------------------------------
void
LAPACK_IMPL(zgttrf)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *DL,
                    DOUBLE_COMPLEX   *D,
                    DOUBLE_COMPLEX   *DU,
                    DOUBLE_COMPLEX   *DU2,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zgttrs --------------------------------------------------------------------
void
LAPACK_IMPL(zgttrs)(const char               *TRANS,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *DL,
                    const DOUBLE_COMPLEX     *D,
                    const DOUBLE_COMPLEX     *DU,
                    const DOUBLE_COMPLEX     *DU2,
                    const INTEGER            *IPIV,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zgtts2 --------------------------------------------------------------------
void
LAPACK_IMPL(zgtts2)(const INTEGER            *ITRANS,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *DL,
                    const DOUBLE_COMPLEX     *D,
                    const DOUBLE_COMPLEX     *DU,
                    const DOUBLE_COMPLEX     *DU2,
                    const INTEGER            *IPIV,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB);

//-- zhbev ---------------------------------------------------------------------
void
LAPACK_IMPL(zhbev)(const char           *JOBZ,
                   const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *KD,
                   DOUBLE_COMPLEX       *AB,
                   const INTEGER        *LDAB,
                   DOUBLE               *W,
                   DOUBLE_COMPLEX       *Z,
                   const INTEGER        *LDZ,
                   DOUBLE_COMPLEX       *WORK,
                   DOUBLE               *RWORK,
                   INTEGER              *INFO);

//-- zhbevd --------------------------------------------------------------------
void
LAPACK_IMPL(zhbevd)(const char       *JOBZ,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *W,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    const INTEGER    *LRWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- zhbevx --------------------------------------------------------------------
void
LAPACK_IMPL(zhbevx)(const char       *JOBZ,
                    const char       *RANGE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *IWORK,
                    INTEGER          *IFAIL,
                    INTEGER          *INFO);

//-- zhbgst --------------------------------------------------------------------
void
LAPACK_IMPL(zhbgst)(const char               *VECT,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *KA,
                    const INTEGER            *KB,
                    DOUBLE_COMPLEX           *AB,
                    const INTEGER            *LDAB,
                    const DOUBLE_COMPLEX     *BB,
                    const INTEGER            *LDBB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zhbgv ---------------------------------------------------------------------
void
LAPACK_IMPL(zhbgv)(const char           *JOBZ,
                   const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *KA,
                   const INTEGER        *KB,
                   DOUBLE_COMPLEX       *AB,
                   const INTEGER        *LDAB,
                   DOUBLE_COMPLEX       *BB,
                   const INTEGER        *LDBB,
                   DOUBLE               *W,
                   DOUBLE_COMPLEX       *Z,
                   const INTEGER        *LDZ,
                   DOUBLE_COMPLEX       *WORK,
                   DOUBLE               *RWORK,
                   INTEGER              *INFO);

//-- zhbgvd --------------------------------------------------------------------
void
LAPACK_IMPL(zhbgvd)(const char       *JOBZ,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KA,
                    const INTEGER    *KB,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    DOUBLE_COMPLEX   *BB,
                    const INTEGER    *LDBB,
                    DOUBLE           *W,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    const INTEGER    *LRWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- zhbgvx --------------------------------------------------------------------
void
LAPACK_IMPL(zhbgvx)(const char       *JOBZ,
                    const char       *RANGE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KA,
                    const INTEGER    *KB,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    DOUBLE_COMPLEX   *BB,
                    const INTEGER    *LDBB,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *IWORK,
                    INTEGER          *IFAIL,
                    INTEGER          *INFO);

//-- zhbtrd --------------------------------------------------------------------
void
LAPACK_IMPL(zhbtrd)(const char       *VECT,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zhecon --------------------------------------------------------------------
void
LAPACK_IMPL(zhecon)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const INTEGER            *IPIV,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zheequb -------------------------------------------------------------------
void
LAPACK_IMPL(zheequb)(const char               *UPLO,
                     const INTEGER            *N,
                     const DOUBLE_COMPLEX     *A,
                     const INTEGER            *LDA,
                     DOUBLE                   *S,
                     DOUBLE                   *SCOND,
                     DOUBLE                   *AMAX,
                     const DOUBLE_COMPLEX     *WORK,
                     INTEGER                  *INFO);

//-- zheev ---------------------------------------------------------------------
void
LAPACK_IMPL(zheev)(const char           *JOBZ,
                   const char           *UPLO,
                   const INTEGER        *N,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   DOUBLE               *W,
                   DOUBLE_COMPLEX       *WORK,
                   const INTEGER        *LWORK,
                   DOUBLE               *RWORK,
                   INTEGER              *INFO);

//-- zheevd --------------------------------------------------------------------
void
LAPACK_IMPL(zheevd)(const char       *JOBZ,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE           *W,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    const INTEGER    *LRWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- zheevr --------------------------------------------------------------------
void
LAPACK_IMPL(zheevr)(const char       *JOBZ,
                    const char       *RANGE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *ISUPPZ,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    const INTEGER    *LRWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- zheevx --------------------------------------------------------------------
void
LAPACK_IMPL(zheevx)(const char       *JOBZ,
                    const char       *RANGE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    INTEGER          *IWORK,
                    INTEGER          *IFAIL,
                    INTEGER          *INFO);

//-- zhegs2 --------------------------------------------------------------------
void
LAPACK_IMPL(zhegs2)(const INTEGER            *ITYPE,
                    const char               *UPLO,
                    const INTEGER            *N,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zhegst --------------------------------------------------------------------
void
LAPACK_IMPL(zhegst)(const INTEGER            *ITYPE,
                    const char               *UPLO,
                    const INTEGER            *N,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zhegv ---------------------------------------------------------------------
void
LAPACK_IMPL(zhegv)(const INTEGER        *ITYPE,
                   const char           *JOBZ,
                   const char           *UPLO,
                   const INTEGER        *N,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   DOUBLE               *W,
                   DOUBLE_COMPLEX       *WORK,
                   const INTEGER        *LWORK,
                   DOUBLE               *RWORK,
                   INTEGER              *INFO);

//-- zhegvd --------------------------------------------------------------------
void
LAPACK_IMPL(zhegvd)(const INTEGER    *ITYPE,
                    const char       *JOBZ,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE           *W,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    const INTEGER    *LRWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- zhegvx --------------------------------------------------------------------
void
LAPACK_IMPL(zhegvx)(const INTEGER    *ITYPE,
                    const char       *JOBZ,
                    const char       *RANGE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    INTEGER          *IWORK,
                    INTEGER          *IFAIL,
                    INTEGER          *INFO);

//-- zherfs --------------------------------------------------------------------
void
LAPACK_IMPL(zherfs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *AF,
                    const INTEGER            *LDAF,
                    const INTEGER            *IPIV,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zherfsx -------------------------------------------------------------------
/*
void
LAPACK_IMPL(zherfsx)(const char               *UPLO,
                     const char               *EQUED,
                     const INTEGER            *N,
                     const INTEGER            *NRHS,
                     const DOUBLE_COMPLEX     *A,
                     const INTEGER            *LDA,
                     const DOUBLE_COMPLEX     *AF,
                     const INTEGER            *LDAF,
                     const INTEGER            *IPIV,
                     DOUBLE                   *S,
                     const DOUBLE_COMPLEX     *B,
                     const INTEGER            *LDB,
                     DOUBLE_COMPLEX           *X,
                     const INTEGER            *LDX,
                     DOUBLE                   *RCOND,
                     DOUBLE                   *BERR,
                     const INTEGER            *N_ERR_BNDS,
                     DOUBLE                   *ERR_BNDS_NORM,
                     DOUBLE                   *ERR_BNDS_COMP,
                     const INTEGER            *NPARAMS,
                     DOUBLE                   *PARAMS,
                     DOUBLE_COMPLEX           *WORK,
                     DOUBLE                   *RWORK,
                     INTEGER                  *INFO);
*/

//-- zhesv ---------------------------------------------------------------------
void
LAPACK_IMPL(zhesv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   INTEGER              *IPIV,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   DOUBLE_COMPLEX       *WORK,
                   const INTEGER        *LWORK,
                   INTEGER              *INFO);

//-- zhesvx --------------------------------------------------------------------
void
LAPACK_IMPL(zhesvx)(const char               *FACT,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *AF,
                    const INTEGER            *LDAF,
                    INTEGER                  *IPIV,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *RCOND,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zhesvxx -------------------------------------------------------------------
/*
void
LAPACK_IMPL(zhesvxx)(const char       *FACT,
                     const char       *UPLO,
                     const INTEGER    *N,
                     const INTEGER    *NRHS,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     DOUBLE_COMPLEX   *AF,
                     const INTEGER    *LDAF,
                     INTEGER          *IPIV,
                     char             *EQUED,
                     DOUBLE           *S,
                     DOUBLE_COMPLEX   *B,
                     const INTEGER    *LDB,
                     DOUBLE_COMPLEX   *X,
                     const INTEGER    *LDX,
                     DOUBLE           *RCOND,
                     DOUBLE           *RPVGRW,
                     DOUBLE           *BERR,
                     const INTEGER    *N_ERR_BNDS,
                     DOUBLE           *ERR_BNDS_NORM,
                     DOUBLE           *ERR_BNDS_COMP,
                     const INTEGER    *NPARAMS,
                     DOUBLE           *PARAMS,
                     DOUBLE_COMPLEX   *WORK,
                     DOUBLE           *RWORK,
                     INTEGER          *INFO);
*/
//-- zheswapr ------------------------------------------------------------------
void
LAPACK_IMPL(zheswapr)(const char           *UPLO,
                      const INTEGER        *N,
                      DOUBLE_COMPLEX       *A,
                      const INTEGER        *LDA,
                      const INTEGER        *I1,
                      const INTEGER        *I2);

//-- zhetd2 --------------------------------------------------------------------
void
LAPACK_IMPL(zhetd2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *TAU,
                    INTEGER          *INFO);

//-- zhetf2 --------------------------------------------------------------------
void
LAPACK_IMPL(zhetf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zhetrd --------------------------------------------------------------------
void
LAPACK_IMPL(zhetrd)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zhetrf --------------------------------------------------------------------
void
LAPACK_IMPL(zhetrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zhetri --------------------------------------------------------------------
void
LAPACK_IMPL(zhetri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zhetri2 -------------------------------------------------------------------
void
LAPACK_IMPL(zhetri2)(const char       *UPLO,
                     const INTEGER    *N,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE_COMPLEX   *WORK,
                     const INTEGER    *LWORK,
                     INTEGER          *INFO);

//-- zhetri2x ------------------------------------------------------------------
void
LAPACK_IMPL(zhetri2x)(const char           *UPLO,
                      const INTEGER        *N,
                      DOUBLE_COMPLEX       *A,
                      const INTEGER        *LDA,
                      const INTEGER        *IPIV,
                      DOUBLE_COMPLEX       *WORK,
                      const INTEGER        *NB,
                      INTEGER              *INFO);

//-- zhetrs --------------------------------------------------------------------
void
LAPACK_IMPL(zhetrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const INTEGER            *IPIV,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zhetrs2 -------------------------------------------------------------------
void
LAPACK_IMPL(zhetrs2)(const char       *UPLO,
                     const INTEGER    *N,
                     const INTEGER    *NRHS,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE_COMPLEX   *B,
                     const INTEGER    *LDB,
                     DOUBLE_COMPLEX   *WORK,
                     INTEGER          *INFO);

//-- zhfrk ---------------------------------------------------------------------
void
LAPACK_IMPL(zhfrk)(const char               *TRANSR,
                   const char               *UPLO,
                   const char               *TRANS,
                   const INTEGER            *N,
                   const INTEGER            *K,
                   const DOUBLE             *ALPHA,
                   const DOUBLE_COMPLEX     *A,
                   const INTEGER            *LDA,
                   const DOUBLE             *BETA,
                   DOUBLE_COMPLEX           *C);

//-- zhgeqz --------------------------------------------------------------------
void
LAPACK_IMPL(zhgeqz)(const char       *JOB,
                    const char       *COMPQ,
                    const char       *COMPZ,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE_COMPLEX   *H,
                    const INTEGER    *LDH,
                    DOUBLE_COMPLEX   *T,
                    const INTEGER    *LDT,
                    DOUBLE_COMPLEX   *ALPHA,
                    DOUBLE_COMPLEX   *BETA,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zhpcon --------------------------------------------------------------------
void
LAPACK_IMPL(zhpcon)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    const INTEGER            *IPIV,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zhpev ---------------------------------------------------------------------
void
LAPACK_IMPL(zhpev)(const char           *JOBZ,
                   const char           *UPLO,
                   const INTEGER        *N,
                   DOUBLE_COMPLEX       *AP,
                   DOUBLE               *W,
                   DOUBLE_COMPLEX       *Z,
                   const INTEGER        *LDZ,
                   DOUBLE_COMPLEX       *WORK,
                   DOUBLE               *RWORK,
                   INTEGER              *INFO);

//-- zhpevd --------------------------------------------------------------------
void
LAPACK_IMPL(zhpevd)(const char       *JOBZ,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    DOUBLE           *W,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    const INTEGER    *LRWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- zhpevx --------------------------------------------------------------------
void
LAPACK_IMPL(zhpevx)(const char       *JOBZ,
                    const char       *RANGE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *IWORK,
                    INTEGER          *IFAIL,
                    INTEGER          *INFO);

//-- zhpgst --------------------------------------------------------------------
void
LAPACK_IMPL(zhpgst)(const INTEGER            *ITYPE,
                    const char               *UPLO,
                    const INTEGER            *N,
                    DOUBLE_COMPLEX           *AP,
                    const DOUBLE_COMPLEX     *BP,
                    INTEGER                  *INFO);

//-- zhpgv ---------------------------------------------------------------------
void
LAPACK_IMPL(zhpgv)(const INTEGER        *ITYPE,
                   const char           *JOBZ,
                   const char           *UPLO,
                   const INTEGER        *N,
                   DOUBLE_COMPLEX       *AP,
                   DOUBLE_COMPLEX       *BP,
                   DOUBLE               *W,
                   DOUBLE_COMPLEX       *Z,
                   const INTEGER        *LDZ,
                   DOUBLE_COMPLEX       *WORK,
                   DOUBLE               *RWORK,
                   INTEGER              *INFO);

//-- zhpgvd --------------------------------------------------------------------
void
LAPACK_IMPL(zhpgvd)(const INTEGER    *ITYPE,
                    const char       *JOBZ,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    DOUBLE_COMPLEX   *BP,
                    DOUBLE           *W,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    const INTEGER    *LRWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- zhpgvx --------------------------------------------------------------------
void
LAPACK_IMPL(zhpgvx)(const INTEGER    *ITYPE,
                    const char       *JOBZ,
                    const char       *RANGE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    DOUBLE_COMPLEX   *BP,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *IWORK,
                    INTEGER          *IFAIL,
                    INTEGER          *INFO);

//-- zhprfs --------------------------------------------------------------------
void
LAPACK_IMPL(zhprfs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    const DOUBLE_COMPLEX     *AFP,
                    const INTEGER            *IPIV,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zhpsv ---------------------------------------------------------------------
void
LAPACK_IMPL(zhpsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *AP,
                   INTEGER              *IPIV,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zhpsvx --------------------------------------------------------------------
void
LAPACK_IMPL(zhpsvx)(const char               *FACT,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE_COMPLEX           *AFP,
                    INTEGER                  *IPIV,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *RCOND,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zhptrd --------------------------------------------------------------------
void
LAPACK_IMPL(zhptrd)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *TAU,
                    INTEGER          *INFO);

//-- zhptrf --------------------------------------------------------------------
void
LAPACK_IMPL(zhptrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zhptri --------------------------------------------------------------------
void
LAPACK_IMPL(zhptri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    const INTEGER    *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zhptrs --------------------------------------------------------------------
void
LAPACK_IMPL(zhptrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    const INTEGER            *IPIV,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zhsein --------------------------------------------------------------------
void
LAPACK_IMPL(zhsein)(const char               *SIDE,
                    const char               *EIGSRC,
                    const char               *INITV,
                    const LOGICAL            *SELECT,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *H,
                    const INTEGER            *LDH,
                    DOUBLE_COMPLEX           *W,
                    DOUBLE_COMPLEX           *VL,
                    const INTEGER            *LDVL,
                    DOUBLE_COMPLEX           *VR,
                    const INTEGER            *LDVR,
                    const INTEGER            *MM,
                    INTEGER                  *M,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *IFAILL,
                    INTEGER                  *IFAILR,
                    INTEGER                  *INFO);

//-- zhseqr --------------------------------------------------------------------
void
LAPACK_IMPL(zhseqr)(const char       *JOB,
                    const char       *COMPZ,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE_COMPLEX   *H,
                    const INTEGER    *LDH,
                    DOUBLE_COMPLEX   *W,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zla_gbamv -----------------------------------------------------------------
void
LAPACK_IMPL(zla_gbamv)(const INTEGER            *TRANS,
                       const INTEGER            *M,
                       const INTEGER            *N,
                       const INTEGER            *KL,
                       const INTEGER            *KU,
                       const DOUBLE             *ALPHA,
                       const DOUBLE_COMPLEX     *AB,
                       const INTEGER            *LDAB,
                       const DOUBLE_COMPLEX     *X,
                       const INTEGER            *INCX,
                       const DOUBLE             *BETA,
                       DOUBLE                   *Y,
                       const INTEGER            *INCY);

//-- zla_gbrcond_c -------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zla_gbrcond_c)(const char               *TRANS,
                           const INTEGER            *N,
                           const INTEGER            *KL,
                           const INTEGER            *KU,
                           const DOUBLE_COMPLEX     *AB,
                           const INTEGER            *LDAB,
                           const DOUBLE_COMPLEX     *AFB,
                           const INTEGER            *LDAFB,
                           const INTEGER            *IPIV,
                           const DOUBLE             *C,
                           const LOGICAL            *CAPPLY,
                           INTEGER                  *INFO,
                           const DOUBLE_COMPLEX     *WORK,
                           const DOUBLE             *RWORK);

//-- zla_gbrcond_x -------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zla_gbrcond_x)(const char               *TRANS,
                           const INTEGER            *N,
                           const INTEGER            *KL,
                           const INTEGER            *KU,
                           const DOUBLE_COMPLEX     *AB,
                           const INTEGER            *LDAB,
                           const DOUBLE_COMPLEX     *AFB,
                           const INTEGER            *LDAFB,
                           const INTEGER            *IPIV,
                           const DOUBLE_COMPLEX     *X,
                           INTEGER                  *INFO,
                           const DOUBLE_COMPLEX     *WORK,
                           const DOUBLE             *RWORK);

//-- zla_gbrfsx_extended -------------------------------------------------------
/*
void
LAPACK_IMPL(zla_gbrfsx_extended)(const INTEGER            *PREC_TYPE,
                                 const INTEGER            *TRANS_TYPE,
                                 const INTEGER            *N,
                                 const INTEGER            *KL,
                                 const INTEGER            *KU,
                                 const INTEGER            *NRHS,
                                 const DOUBLE_COMPLEX     *AB,
                                 const INTEGER            *LDAB,
                                 const DOUBLE_COMPLEX     *AFB,
                                 const INTEGER            *LDAFB,
                                 const INTEGER            *IPIV,
                                 const LOGICAL            *COLEQU,
                                 const DOUBLE             *C,
                                 const DOUBLE_COMPLEX     *B,
                                 const INTEGER            *LDB,
                                 DOUBLE_COMPLEX           *Y,
                                 const INTEGER            *LDY,
                                 DOUBLE                   *BERR_OUT,
                                 const INTEGER            *N_NORMS,
                                 DOUBLE                   *ERR_BNDS_NORM,
                                 DOUBLE                   *ERR_BNDS_COMP,
                                 const DOUBLE_COMPLEX     *RES,
                                 const DOUBLE             *AYB,
                                 const DOUBLE_COMPLEX     *DY,
                                 const DOUBLE_COMPLEX     *Y_TAIL,
                                 const DOUBLE             *RCOND,
                                 const INTEGER            *ITHRESH,
                                 const DOUBLE             *RTHRESH,
                                 const DOUBLE             *DZ_UB,
                                 const LOGICAL            *IGNORE_CWISE,
                                 INTEGER                  *INFO);
*/

//-- zla_gbrpvgrw --------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zla_gbrpvgrw)(const INTEGER            *N,
                          const INTEGER            *KL,
                          const INTEGER            *KU,
                          const INTEGER            *NCOLS,
                          const DOUBLE_COMPLEX     *AB,
                          const INTEGER            *LDAB,
                          const DOUBLE_COMPLEX     *AFB,
                          const INTEGER            *LDAFB);

//-- zla_geamv -----------------------------------------------------------------
void
LAPACK_IMPL(zla_geamv)(const INTEGER            *TRANS,
                       const INTEGER            *M,
                       const INTEGER            *N,
                       const DOUBLE             *ALPHA,
                       const DOUBLE_COMPLEX     *A,
                       const INTEGER            *LDA,
                       const DOUBLE_COMPLEX     *X,
                       const INTEGER            *INCX,
                       const DOUBLE             *BETA,
                       DOUBLE                   *Y,
                       const INTEGER            *INCY);

//-- zla_gercond_c -------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zla_gercond_c)(const char               *TRANS,
                           const INTEGER            *N,
                           const DOUBLE_COMPLEX     *A,
                           const INTEGER            *LDA,
                           const DOUBLE_COMPLEX     *AF,
                           const INTEGER            *LDAF,
                           const INTEGER            *IPIV,
                           const DOUBLE             *C,
                           const LOGICAL            *CAPPLY,
                           INTEGER                  *INFO,
                           const DOUBLE_COMPLEX     *WORK,
                           const DOUBLE             *RWORK);

//-- zla_gercond_x -------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zla_gercond_x)(const char               *TRANS,
                           const INTEGER            *N,
                           const DOUBLE_COMPLEX     *A,
                           const INTEGER            *LDA,
                           const DOUBLE_COMPLEX     *AF,
                           const INTEGER            *LDAF,
                           const INTEGER            *IPIV,
                           const DOUBLE_COMPLEX     *X,
                           INTEGER                  *INFO,
                           const DOUBLE_COMPLEX     *WORK,
                           const DOUBLE             *RWORK);

//-- zla_gerfsx_extended -------------------------------------------------------
/*
void
LAPACK_IMPL(zla_gerfsx_extended)(const INTEGER            *PREC_TYPE,
                                 const INTEGER            *TRANS_TYPE,
                                 const INTEGER            *N,
                                 const INTEGER            *NRHS,
                                 const DOUBLE_COMPLEX     *A,
                                 const INTEGER            *LDA,
                                 const DOUBLE_COMPLEX     *AF,
                                 const INTEGER            *LDAF,
                                 const INTEGER            *IPIV,
                                 const LOGICAL            *COLEQU,
                                 const DOUBLE             *C,
                                 const DOUBLE_COMPLEX     *B,
                                 const INTEGER            *LDB,
                                 DOUBLE_COMPLEX           *Y,
                                 const INTEGER            *LDY,
                                 DOUBLE                   *BERR_OUT,
                                 const INTEGER            *N_NORMS,
                                 const DOUBLE             *ERRS_N,
                                 const DOUBLE             *ERRS_C,
                                 const DOUBLE_COMPLEX     *RES,
                                 const DOUBLE             *AYB,
                                 const DOUBLE_COMPLEX     *DY,
                                 const DOUBLE_COMPLEX     *Y_TAIL,
                                 const DOUBLE             *RCOND,
                                 const INTEGER            *ITHRESH,
                                 const DOUBLE             *RTHRESH,
                                 const DOUBLE             *DZ_UB,
                                 const LOGICAL            *IGNORE_CWISE,
                                 INTEGER                  *INFO);
*/

//-- zla_heamv -----------------------------------------------------------------
void
LAPACK_IMPL(zla_heamv)(const INTEGER            *UPLO,
                       const INTEGER            *N,
                       const DOUBLE             *ALPHA,
                       const DOUBLE_COMPLEX     *A,
                       const INTEGER            *LDA,
                       const DOUBLE_COMPLEX     *X,
                       const INTEGER            *INCX,
                       const DOUBLE             *BETA,
                       DOUBLE                   *Y,
                       const INTEGER            *INCY);

//-- zla_hercond_c -------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zla_hercond_c)(const char               *UPLO,
                           const INTEGER            *N,
                           const DOUBLE_COMPLEX     *A,
                           const INTEGER            *LDA,
                           const DOUBLE_COMPLEX     *AF,
                           const INTEGER            *LDAF,
                           const INTEGER            *IPIV,
                           const DOUBLE             *C,
                           const LOGICAL            *CAPPLY,
                           INTEGER                  *INFO,
                           const DOUBLE_COMPLEX     *WORK,
                           const DOUBLE             *RWORK);

//-- zla_hercond_x -------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zla_hercond_x)(const char               *UPLO,
                           const INTEGER            *N,
                           const DOUBLE_COMPLEX     *A,
                           const INTEGER            *LDA,
                           const DOUBLE_COMPLEX     *AF,
                           const INTEGER            *LDAF,
                           const INTEGER            *IPIV,
                           const DOUBLE_COMPLEX     *X,
                           INTEGER                  *INFO,
                           const DOUBLE_COMPLEX     *WORK,
                           const DOUBLE             *RWORK);

//-- zla_herfsx_extended -------------------------------------------------------
/*
void
LAPACK_IMPL(zla_herfsx_extended)(const INTEGER            *PREC_TYPE,
                                 const char               *UPLO,
                                 const INTEGER            *N,
                                 const INTEGER            *NRHS,
                                 const DOUBLE_COMPLEX     *A,
                                 const INTEGER            *LDA,
                                 const DOUBLE_COMPLEX     *AF,
                                 const INTEGER            *LDAF,
                                 const INTEGER            *IPIV,
                                 const LOGICAL            *COLEQU,
                                 const DOUBLE             *C,
                                 const DOUBLE_COMPLEX     *B,
                                 const INTEGER            *LDB,
                                 DOUBLE_COMPLEX           *Y,
                                 const INTEGER            *LDY,
                                 DOUBLE                   *BERR_OUT,
                                 const INTEGER            *N_NORMS,
                                 DOUBLE                   *ERR_BNDS_NORM,
                                 DOUBLE                   *ERR_BNDS_COMP,
                                 const DOUBLE_COMPLEX     *RES,
                                 const DOUBLE             *AYB,
                                 const DOUBLE_COMPLEX     *DY,
                                 const DOUBLE_COMPLEX     *Y_TAIL,
                                 const DOUBLE             *RCOND,
                                 const INTEGER            *ITHRESH,
                                 const DOUBLE             *RTHRESH,
                                 const DOUBLE             *DZ_UB,
                                 const LOGICAL            *IGNORE_CWISE,
                                 INTEGER                  *INFO);
*/

//-- zla_herpvgrw --------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zla_herpvgrw)(const char               *UPLO,
                          const INTEGER            *N,
                          const INTEGER            *INFO,
                          const DOUBLE_COMPLEX     *A,
                          const INTEGER            *LDA,
                          const DOUBLE_COMPLEX     *AF,
                          const INTEGER            *LDAF,
                          const INTEGER            *IPIV,
                          const DOUBLE             *WORK);

//-- zla_lin_berr --------------------------------------------------------------
void
LAPACK_IMPL(zla_lin_berr)(const INTEGER            *N,
                          const INTEGER            *NZ,
                          const INTEGER            *NRHS,
                          const DOUBLE_COMPLEX     *RES,
                          const DOUBLE             *AYB,
                          DOUBLE                   *BERR);

//-- zla_porcond_c -------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zla_porcond_c)(const char               *UPLO,
                           const INTEGER            *N,
                           const DOUBLE_COMPLEX     *A,
                           const INTEGER            *LDA,
                           const DOUBLE_COMPLEX     *AF,
                           const INTEGER            *LDAF,
                           const DOUBLE             *C,
                           const LOGICAL            *CAPPLY,
                           INTEGER                  *INFO,
                           const DOUBLE_COMPLEX     *WORK,
                           const DOUBLE             *RWORK);

//-- zla_porcond_x -------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zla_porcond_x)(const char               *UPLO,
                           const INTEGER            *N,
                           const DOUBLE_COMPLEX     *A,
                           const INTEGER            *LDA,
                           const DOUBLE_COMPLEX     *AF,
                           const INTEGER            *LDAF,
                           const DOUBLE_COMPLEX     *X,
                           INTEGER                  *INFO,
                           const DOUBLE_COMPLEX     *WORK,
                           const DOUBLE             *RWORK);

//-- zla_porfsx_extended -------------------------------------------------------
/*
void
LAPACK_IMPL(zla_porfsx_extended)(const INTEGER            *PREC_TYPE,
                                 const char               *UPLO,
                                 const INTEGER            *N,
                                 const INTEGER            *NRHS,
                                 const DOUBLE_COMPLEX     *A,
                                 const INTEGER            *LDA,
                                 const DOUBLE_COMPLEX     *AF,
                                 const INTEGER            *LDAF,
                                 const LOGICAL            *COLEQU,
                                 const DOUBLE             *C,
                                 const DOUBLE_COMPLEX     *B,
                                 const INTEGER            *LDB,
                                 DOUBLE_COMPLEX           *Y,
                                 const INTEGER            *LDY,
                                 DOUBLE                   *BERR_OUT,
                                 const INTEGER            *N_NORMS,
                                 DOUBLE                   *ERR_BNDS_NORM,
                                 DOUBLE                   *ERR_BNDS_COMP,
                                 const DOUBLE_COMPLEX     *RES,
                                 const DOUBLE             *AYB,
                                 const DOUBLE_COMPLEX     *DY,
                                 const DOUBLE_COMPLEX     *Y_TAIL,
                                 const DOUBLE             *RCOND,
                                 const INTEGER            *ITHRESH,
                                 const DOUBLE             *RTHRESH,
                                 const DOUBLE             *DZ_UB,
                                 const LOGICAL            *IGNORE_CWISE,
                                 INTEGER                  *INFO);
*/

//-- zla_porpvgrw --------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zla_porpvgrw)(const char               *UPLO,
                          const INTEGER            *NCOLS,
                          const DOUBLE_COMPLEX     *A,
                          const INTEGER            *LDA,
                          const DOUBLE_COMPLEX     *AF,
                          const INTEGER            *LDAF,
                          const DOUBLE             *WORK);

//-- zla_rpvgrw ----------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zla_rpvgrw)(const INTEGER            *N,
                        const INTEGER            *NCOLS,
                        const DOUBLE_COMPLEX     *A,
                        const INTEGER            *LDA,
                        const DOUBLE_COMPLEX     *AF,
                        const INTEGER            *LDAF);

//-- zla_syamv -----------------------------------------------------------------
void
LAPACK_IMPL(zla_syamv)(const INTEGER            *UPLO,
                       const INTEGER            *N,
                       const DOUBLE             *ALPHA,
                       const DOUBLE_COMPLEX     *A,
                       const INTEGER            *LDA,
                       const DOUBLE_COMPLEX     *X,
                       const INTEGER            *INCX,
                       const DOUBLE             *BETA,
                       DOUBLE                   *Y,
                       const INTEGER            *INCY);

//-- zla_syrcond_c -------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zla_syrcond_c)(const char               *UPLO,
                           const INTEGER            *N,
                           const DOUBLE_COMPLEX     *A,
                           const INTEGER            *LDA,
                           const DOUBLE_COMPLEX     *AF,
                           const INTEGER            *LDAF,
                           const INTEGER            *IPIV,
                           const DOUBLE             *C,
                           const LOGICAL            *CAPPLY,
                           INTEGER                  *INFO,
                           const DOUBLE_COMPLEX     *WORK,
                           const DOUBLE             *RWORK);

//-- zla_syrcond_x -------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zla_syrcond_x)(const char               *UPLO,
                           const INTEGER            *N,
                           const DOUBLE_COMPLEX     *A,
                           const INTEGER            *LDA,
                           const DOUBLE_COMPLEX     *AF,
                           const INTEGER            *LDAF,
                           const INTEGER            *IPIV,
                           const DOUBLE_COMPLEX     *X,
                           INTEGER                  *INFO,
                           const DOUBLE_COMPLEX     *WORK,
                           const DOUBLE             *RWORK);

//-- zla_syrfsx_extended -------------------------------------------------------
/*
void
LAPACK_IMPL(zla_syrfsx_extended)(const INTEGER            *PREC_TYPE,
                                 const char               *UPLO,
                                 const INTEGER            *N,
                                 const INTEGER            *NRHS,
                                 const DOUBLE_COMPLEX     *A,
                                 const INTEGER            *LDA,
                                 const DOUBLE_COMPLEX     *AF,
                                 const INTEGER            *LDAF,
                                 const INTEGER            *IPIV,
                                 const LOGICAL            *COLEQU,
                                 const DOUBLE             *C,
                                 const DOUBLE_COMPLEX     *B,
                                 const INTEGER            *LDB,
                                 DOUBLE_COMPLEX           *Y,
                                 const INTEGER            *LDY,
                                 DOUBLE                   *BERR_OUT,
                                 const INTEGER            *N_NORMS,
                                 DOUBLE                   *ERR_BNDS_NORM,
                                 DOUBLE                   *ERR_BNDS_COMP,
                                 const DOUBLE_COMPLEX     *RES,
                                 const DOUBLE             *AYB,
                                 const DOUBLE_COMPLEX     *DY,
                                 const DOUBLE_COMPLEX     *Y_TAIL,
                                 const DOUBLE             *RCOND,
                                 const INTEGER            *ITHRESH,
                                 const DOUBLE             *RTHRESH,
                                 const DOUBLE             *DZ_UB,
                                 const LOGICAL            *IGNORE_CWISE,
                                 INTEGER                  *INFO);
*/

//-- zla_syrpvgrw --------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zla_syrpvgrw)(const char               *UPLO,
                          const INTEGER            *N,
                          const INTEGER            *INFO,
                          const DOUBLE_COMPLEX     *A,
                          const INTEGER            *LDA,
                          const DOUBLE_COMPLEX     *AF,
                          const INTEGER            *LDAF,
                          const INTEGER            *IPIV,
                          const DOUBLE             *WORK);

//-- zla_wwaddw ----------------------------------------------------------------
void
LAPACK_IMPL(zla_wwaddw)(const INTEGER            *N,
                        DOUBLE_COMPLEX           *X,
                        DOUBLE_COMPLEX           *Y,
                        const DOUBLE_COMPLEX     *W);

//-- zlabrd --------------------------------------------------------------------
void
LAPACK_IMPL(zlabrd)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *NB,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *TAUQ,
                    DOUBLE_COMPLEX   *TAUP,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *LDX,
                    DOUBLE_COMPLEX   *Y,
                    const INTEGER    *LDY);

//-- zlacgv --------------------------------------------------------------------
void
LAPACK_IMPL(zlacgv)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *INCX);

//-- zlacn2 --------------------------------------------------------------------
void
LAPACK_IMPL(zlacn2)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *V,
                    DOUBLE_COMPLEX   *X,
                    DOUBLE           *EST,
                    INTEGER          *KASE,
                    INTEGER          *ISAVE);

//-- zlacon --------------------------------------------------------------------
void
LAPACK_IMPL(zlacon)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *V,
                    DOUBLE_COMPLEX   *X,
                    DOUBLE           *EST,
                    INTEGER          *KASE);

//-- zlacp2 --------------------------------------------------------------------
void
LAPACK_IMPL(zlacp2)(const char       *UPLO,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB);

//-- zlacpy --------------------------------------------------------------------
void
LAPACK_IMPL(zlacpy)(const char               *UPLO,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB);

//-- zlacrm --------------------------------------------------------------------
void
LAPACK_IMPL(zlacrm)(const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE             *B,
                    const INTEGER            *LDB,
                    const DOUBLE_COMPLEX     *C,
                    const INTEGER            *LDC,
                    DOUBLE                   *RWORK);

//-- zlacrt --------------------------------------------------------------------
void
LAPACK_IMPL(zlacrt)(const INTEGER            *N,
                    DOUBLE_COMPLEX           *CX,
                    const INTEGER            *INCX,
                    DOUBLE_COMPLEX           *CY,
                    const INTEGER            *INCY,
                    const DOUBLE_COMPLEX     *C,
                    const DOUBLE_COMPLEX     *S);

//-- zladiv --------------------------------------------------------------------
/*
void
LAPACK_IMPL(zladiv)(DOUBLE_COMPLEX           *ret,
                    const DOUBLE_COMPLEX     *X,
                    const DOUBLE_COMPLEX     *Y);
*/

//-- zlaed0 --------------------------------------------------------------------
void
LAPACK_IMPL(zlaed0)(const INTEGER    *QSIZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    DOUBLE_COMPLEX   *QSTORE,
                    const INTEGER    *LDQS,
                    DOUBLE           *RWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- zlaed7 --------------------------------------------------------------------
void
LAPACK_IMPL(zlaed7)(const INTEGER    *N,
                    const INTEGER    *CUTPNT,
                    const INTEGER    *QSIZ,
                    const INTEGER    *TLVLS,
                    const INTEGER    *CURLVL,
                    const INTEGER    *CURPBM,
                    DOUBLE           *D,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    const DOUBLE     *RHO,
                    INTEGER          *INDXQ,
                    DOUBLE           *QSTORE,
                    INTEGER          *QPTR,
                    const INTEGER    *PRMPTR,
                    const INTEGER    *PERM,
                    const INTEGER    *GIVPTR,
                    const INTEGER    *GIVCOL,
                    const DOUBLE     *GIVNUM,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- zlaed8 --------------------------------------------------------------------
void
LAPACK_IMPL(zlaed8)(INTEGER          *K,
                    const INTEGER    *N,
                    const INTEGER    *QSIZ,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    DOUBLE           *D,
                    DOUBLE           *RHO,
                    const INTEGER    *CUTPNT,
                    const DOUBLE     *Z,
                    DOUBLE           *DLAMDA,
                    DOUBLE_COMPLEX   *Q2,
                    const INTEGER    *LDQ2,
                    DOUBLE           *W,
                    INTEGER          *INDXP,
                    INTEGER          *INDX,
                    const INTEGER    *INDXQ,
                    INTEGER          *PERM,
                    INTEGER          *GIVPTR,
                    INTEGER          *GIVCOL,
                    DOUBLE           *GIVNUM,
                    INTEGER          *INFO);

//-- zlaein --------------------------------------------------------------------
void
LAPACK_IMPL(zlaein)(const LOGICAL            *RIGHTV,
                    const LOGICAL            *NOINIT,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *H,
                    const INTEGER            *LDH,
                    const DOUBLE_COMPLEX     *W,
                    DOUBLE_COMPLEX           *V,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    DOUBLE                   *RWORK,
                    const DOUBLE             *EPS3,
                    const DOUBLE             *SMLNUM,
                    INTEGER                  *INFO);

//-- zlaesy --------------------------------------------------------------------
void
LAPACK_IMPL(zlaesy)(const DOUBLE_COMPLEX     *A,
                    const DOUBLE_COMPLEX     *B,
                    const DOUBLE_COMPLEX     *C,
                    DOUBLE_COMPLEX           *RT1,
                    DOUBLE_COMPLEX           *RT2,
                    DOUBLE_COMPLEX           *EVSCAL,
                    DOUBLE_COMPLEX           *CS1,
                    DOUBLE_COMPLEX           *SN1);

//-- zlaev2 --------------------------------------------------------------------
void
LAPACK_IMPL(zlaev2)(const DOUBLE_COMPLEX     *A,
                    const DOUBLE_COMPLEX     *B,
                    const DOUBLE_COMPLEX     *C,
                    DOUBLE                   *RT1,
                    DOUBLE                   *RT2,
                    DOUBLE                   *CS1,
                    DOUBLE_COMPLEX           *SN1);

//-- zlag2c --------------------------------------------------------------------
void
LAPACK_IMPL(zlag2c)(const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    FLOAT_COMPLEX            *SA,
                    const INTEGER            *LDSA,
                    INTEGER                  *INFO);

//-- zlags2 --------------------------------------------------------------------
void
LAPACK_IMPL(zlags2)(const LOGICAL            *UPPER,
                    const DOUBLE             *A1,
                    const DOUBLE_COMPLEX     *A2,
                    const DOUBLE             *A3,
                    const DOUBLE             *B1,
                    const DOUBLE_COMPLEX     *B2,
                    const DOUBLE             *B3,
                    DOUBLE                   *CSU,
                    DOUBLE_COMPLEX           *SNU,
                    DOUBLE                   *CSV,
                    DOUBLE_COMPLEX           *SNV,
                    DOUBLE                   *CSQ,
                    DOUBLE_COMPLEX           *SNQ);

//-- zlagtm --------------------------------------------------------------------
void
LAPACK_IMPL(zlagtm)(const char               *TRANS,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE             *ALPHA,
                    const DOUBLE_COMPLEX     *DL,
                    const DOUBLE_COMPLEX     *D,
                    const DOUBLE_COMPLEX     *DU,
                    const DOUBLE_COMPLEX     *X,
                    const INTEGER            *LDX,
                    const DOUBLE             *BETA,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB);

//-- zlahef --------------------------------------------------------------------
void
LAPACK_IMPL(zlahef)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NB,
                    INTEGER          *KB,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    DOUBLE_COMPLEX   *W,
                    const INTEGER    *LDW,
                    INTEGER          *INFO);

//-- zlahqr --------------------------------------------------------------------
void
LAPACK_IMPL(zlahqr)(const LOGICAL    *WANTT,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE_COMPLEX   *H,
                    const INTEGER    *LDH,
                    DOUBLE_COMPLEX   *W,
                    const INTEGER    *ILOZ,
                    const INTEGER    *IHIZ,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *INFO);

//-- zlahr2 --------------------------------------------------------------------
void
LAPACK_IMPL(zlahr2)(const INTEGER    *N,
                    const INTEGER    *K,
                    const INTEGER    *NB,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *T,
                    const INTEGER    *LDT,
                    DOUBLE_COMPLEX   *Y,
                    const INTEGER    *LDY);

//-- zlahrd --------------------------------------------------------------------
void
LAPACK_IMPL(zlahrd)(const INTEGER    *N,
                    const INTEGER    *K,
                    const INTEGER    *NB,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *T,
                    const INTEGER    *LDT,
                    DOUBLE_COMPLEX   *Y,
                    const INTEGER    *LDY);

//-- zlaic1 --------------------------------------------------------------------
void
LAPACK_IMPL(zlaic1)(const INTEGER            *JOB,
                    const INTEGER            *J,
                    const DOUBLE_COMPLEX     *X,
                    const DOUBLE             *SEST,
                    const DOUBLE_COMPLEX     *W,
                    const DOUBLE_COMPLEX     *GAMMA,
                    DOUBLE                   *SESTPR,
                    DOUBLE_COMPLEX           *S,
                    DOUBLE_COMPLEX           *C);

//-- zlals0 --------------------------------------------------------------------
void
LAPACK_IMPL(zlals0)(const INTEGER    *ICOMPQ,
                    const INTEGER    *NL,
                    const INTEGER    *NR,
                    const INTEGER    *SQRE,
                    const INTEGER    *NRHS,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *BX,
                    const INTEGER    *LDBX,
                    const INTEGER    *PERM,
                    const INTEGER    *GIVPTR,
                    const INTEGER    *GIVCOL,
                    const INTEGER    *LDGCOL,
                    const DOUBLE     *GIVNUM,
                    const INTEGER    *LDGNUM,
                    const DOUBLE     *POLES,
                    const DOUBLE     *DIFL,
                    const DOUBLE     *DIFR,
                    const DOUBLE     *Z,
                    const INTEGER    *K,
                    const DOUBLE     *C,
                    const DOUBLE     *S,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zlalsa --------------------------------------------------------------------
void
LAPACK_IMPL(zlalsa)(const INTEGER    *ICOMPQ,
                    const INTEGER    *SMLSIZ,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *BX,
                    const INTEGER    *LDBX,
                    const DOUBLE     *U,
                    const INTEGER    *LDU,
                    const DOUBLE     *VT,
                    const INTEGER    *K,
                    const DOUBLE     *DIFL,
                    const DOUBLE     *DIFR,
                    const DOUBLE     *Z,
                    const DOUBLE     *POLES,
                    const INTEGER    *GIVPTR,
                    const INTEGER    *GIVCOL,
                    const INTEGER    *LDGCOL,
                    const INTEGER    *PERM,
                    const DOUBLE     *GIVNUM,
                    const DOUBLE     *C,
                    const DOUBLE     *S,
                    DOUBLE           *RWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- zlalsd --------------------------------------------------------------------
void
LAPACK_IMPL(zlalsd)(const char       *UPLO,
                    const INTEGER    *SMLSIZ,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    const DOUBLE     *RCOND,
                    INTEGER          *RANK,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- zlangb --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zlangb)(const char               *NORM,
                    const INTEGER            *N,
                    const INTEGER            *KL,
                    const INTEGER            *KU,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE                   *WORK);

//-- zlange --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zlange)(const char               *NORM,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *WORK);

//-- zlangt --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zlangt)(const char               *NORM,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *DL,
                    const DOUBLE_COMPLEX     *D,
                    const DOUBLE_COMPLEX     *DU);

//-- zlanhb --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zlanhb)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE                   *WORK);

//-- zlanhe --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zlanhe)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *WORK);

//-- zlanhf --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zlanhf)(const char               *NORM,
                    const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    DOUBLE                   *WORK);

//-- zlanhp --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zlanhp)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE                   *WORK);

//-- zlanhs --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zlanhs)(const char               *NORM,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *WORK);

//-- zlanht --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zlanht)(const char               *NORM,
                    const INTEGER            *N,
                    const DOUBLE             *D,
                    const DOUBLE_COMPLEX     *E);

//-- zlansb --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zlansb)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE                   *WORK);

//-- zlansp --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zlansp)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE                   *WORK);

//-- zlansy --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zlansy)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *WORK);

//-- zlantb --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zlantb)(const char               *NORM,
                    const char               *UPLO,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE                   *WORK);

//-- zlantp --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zlantp)(const char               *NORM,
                    const char               *UPLO,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE                   *WORK);

//-- zlantr --------------------------------------------------------------------
DOUBLE
LAPACK_IMPL(zlantr)(const char               *NORM,
                    const char               *UPLO,
                    const char               *DIAG,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *WORK);

//-- zlapll --------------------------------------------------------------------
void
LAPACK_IMPL(zlapll)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *INCX,
                    DOUBLE_COMPLEX   *Y,
                    const INTEGER    *INCY,
                    DOUBLE           *SSMIN);

//-- zlapmr --------------------------------------------------------------------
void
LAPACK_IMPL(zlapmr)(const LOGICAL    *FORWRD,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *LDX,
                    INTEGER          *K);

//-- zlapmt --------------------------------------------------------------------
void
LAPACK_IMPL(zlapmt)(const LOGICAL    *FORWRD,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *LDX,
                    INTEGER          *K);

//-- zlaqgb --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqgb)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    const DOUBLE     *R,
                    const DOUBLE     *C,
                    const DOUBLE     *ROWCND,
                    const DOUBLE     *COLCND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- zlaqge --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqge)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *R,
                    const DOUBLE     *C,
                    const DOUBLE     *ROWCND,
                    const DOUBLE     *COLCND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- zlaqhb --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqhb)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- zlaqhe --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqhe)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- zlaqhp --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqhp)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- zlaqp2 --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqp2)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *OFFSET,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *JPVT,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE           *VN1,
                    DOUBLE           *VN2,
                    DOUBLE_COMPLEX   *WORK);

//-- zlaqps --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqps)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *OFFSET,
                    const INTEGER    *NB,
                    INTEGER          *KB,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *JPVT,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE           *VN1,
                    DOUBLE           *VN2,
                    DOUBLE_COMPLEX   *AUXV,
                    DOUBLE_COMPLEX   *F,
                    const INTEGER    *LDF);

//-- zlaqr0 --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqr0)(const LOGICAL    *WANTT,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE_COMPLEX   *H,
                    const INTEGER    *LDH,
                    DOUBLE_COMPLEX   *W,
                    const INTEGER    *ILOZ,
                    const INTEGER    *IHIZ,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zlaqr1 --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqr1)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *H,
                    const INTEGER            *LDH,
                    const DOUBLE_COMPLEX     *S1,
                    const DOUBLE_COMPLEX     *S2,
                    DOUBLE_COMPLEX           *V);

//-- zlaqr2 --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqr2)(const LOGICAL    *WANTT,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *N,
                    const INTEGER    *KTOP,
                    const INTEGER    *KBOT,
                    const INTEGER    *NW,
                    DOUBLE_COMPLEX   *H,
                    const INTEGER    *LDH,
                    const INTEGER    *ILOZ,
                    const INTEGER    *IHIZ,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *NS,
                    INTEGER          *ND,
                    DOUBLE_COMPLEX   *SH,
                    DOUBLE_COMPLEX   *V,
                    const INTEGER    *LDV,
                    const INTEGER    *NH,
                    DOUBLE_COMPLEX   *T,
                    const INTEGER    *LDT,
                    const INTEGER    *NV,
                    DOUBLE_COMPLEX   *WV,
                    const INTEGER    *LDWV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK);

//-- zlaqr3 --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqr3)(const LOGICAL    *WANTT,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *N,
                    const INTEGER    *KTOP,
                    const INTEGER    *KBOT,
                    const INTEGER    *NW,
                    DOUBLE_COMPLEX   *H,
                    const INTEGER    *LDH,
                    const INTEGER    *ILOZ,
                    const INTEGER    *IHIZ,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *NS,
                    INTEGER          *ND,
                    DOUBLE_COMPLEX   *SH,
                    DOUBLE_COMPLEX   *V,
                    const INTEGER    *LDV,
                    const INTEGER    *NH,
                    DOUBLE_COMPLEX   *T,
                    const INTEGER    *LDT,
                    const INTEGER    *NV,
                    DOUBLE_COMPLEX   *WV,
                    const INTEGER    *LDWV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK);

//-- zlaqr4 --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqr4)(const LOGICAL    *WANTT,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE_COMPLEX   *H,
                    const INTEGER    *LDH,
                    DOUBLE_COMPLEX   *W,
                    const INTEGER    *ILOZ,
                    const INTEGER    *IHIZ,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zlaqr5 --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqr5)(const LOGICAL    *WANTT,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *KACC22,
                    const INTEGER    *N,
                    const INTEGER    *KTOP,
                    const INTEGER    *KBOT,
                    const INTEGER    *NSHFTS,
                    DOUBLE_COMPLEX   *S,
                    DOUBLE_COMPLEX   *H,
                    const INTEGER    *LDH,
                    const INTEGER    *ILOZ,
                    const INTEGER    *IHIZ,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *V,
                    const INTEGER    *LDV,
                    DOUBLE_COMPLEX   *U,
                    const INTEGER    *LDU,
                    const INTEGER    *NV,
                    DOUBLE_COMPLEX   *WV,
                    const INTEGER    *LDWV,
                    const INTEGER    *NH,
                    DOUBLE_COMPLEX   *WH,
                    const INTEGER    *LDWH);

//-- zlaqsb --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqsb)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- zlaqsp --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqsp)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- zlaqsy --------------------------------------------------------------------
void
LAPACK_IMPL(zlaqsy)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- zlar1v --------------------------------------------------------------------
void
LAPACK_IMPL(zlar1v)(const INTEGER    *N,
                    const INTEGER    *B1,
                    const INTEGER    *BN,
                    const DOUBLE     *LAMBDA,
                    const DOUBLE     *D,
                    const DOUBLE     *L,
                    const DOUBLE     *LD,
                    const DOUBLE     *LLD,
                    const DOUBLE     *PIVMIN,
                    const DOUBLE     *GAPTOL,
                    DOUBLE_COMPLEX   *Z,
                    const LOGICAL    *WANTNC,
                    INTEGER          *NEGCNT,
                    DOUBLE           *ZTZ,
                    DOUBLE           *MINGMA,
                    INTEGER          *R,
                    INTEGER          *ISUPPZ,
                    DOUBLE           *NRMINV,
                    DOUBLE           *RESID,
                    DOUBLE           *RQCORR,
                    DOUBLE           *WORK);

//-- zlar2v --------------------------------------------------------------------
void
LAPACK_IMPL(zlar2v)(const INTEGER            *N,
                    DOUBLE_COMPLEX           *X,
                    DOUBLE_COMPLEX           *Y,
                    DOUBLE_COMPLEX           *Z,
                    const INTEGER            *INCX,
                    const DOUBLE             *C,
                    const DOUBLE_COMPLEX     *S,
                    const INTEGER            *INCC);

//-- zlarcm --------------------------------------------------------------------
void
LAPACK_IMPL(zlarcm)(const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE             *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    const DOUBLE_COMPLEX     *C,
                    const INTEGER            *LDC,
                    DOUBLE                   *RWORK);

//-- zlarf ---------------------------------------------------------------------
void
LAPACK_IMPL(zlarf)(const char               *SIDE,
                   const INTEGER            *M,
                   const INTEGER            *N,
                   const DOUBLE_COMPLEX     *V,
                   const INTEGER            *INCV,
                   const DOUBLE_COMPLEX     *TAU,
                   DOUBLE_COMPLEX           *C,
                   const INTEGER            *LDC,
                   DOUBLE_COMPLEX           *WORK);

//-- zlarfb --------------------------------------------------------------------
void
LAPACK_IMPL(zlarfb)(const char               *SIDE,
                    const char               *TRANS,
                    const char               *DIRECT,
                    const char               *STOREV,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *V,
                    const INTEGER            *LDV,
                    const DOUBLE_COMPLEX     *T,
                    const INTEGER            *LDT,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LDWORK);

//-- zlarfg --------------------------------------------------------------------
void
LAPACK_IMPL(zlarfg)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *ALPHA,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *INCX,
                    DOUBLE_COMPLEX   *TAU);

//-- zlarfgp -------------------------------------------------------------------
void
LAPACK_IMPL(zlarfgp)(const INTEGER    *N,
                     DOUBLE_COMPLEX   *ALPHA,
                     DOUBLE_COMPLEX   *X,
                     const INTEGER    *INCX,
                     DOUBLE_COMPLEX   *TAU);

//-- zlarft --------------------------------------------------------------------
void
LAPACK_IMPL(zlarft)(const char               *DIRECT,
                    const char               *STOREV,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *V,
                    const INTEGER            *LDV,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *T,
                    const INTEGER            *LDT);

//-- zlarfx --------------------------------------------------------------------
void
LAPACK_IMPL(zlarfx)(const char               *SIDE,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *V,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK);

//-- zlargv --------------------------------------------------------------------
void
LAPACK_IMPL(zlargv)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *INCX,
                    DOUBLE_COMPLEX   *Y,
                    const INTEGER    *INCY,
                    DOUBLE           *C,
                    const INTEGER    *INCC);

//-- zlarnv --------------------------------------------------------------------
void
LAPACK_IMPL(zlarnv)(const INTEGER    *IDIST,
                    INTEGER          *ISEED,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *X);

//-- zlarrv --------------------------------------------------------------------
void
LAPACK_IMPL(zlarrv)(const INTEGER    *N,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    DOUBLE           *D,
                    DOUBLE           *L,
                    const DOUBLE     *PIVMIN,
                    const INTEGER    *ISPLIT,
                    const INTEGER    *M,
                    const INTEGER    *DOL,
                    const INTEGER    *DOU,
                    const DOUBLE     *MINRGP,
                    const DOUBLE     *RTOL1,
                    const DOUBLE     *RTOL2,
                    DOUBLE           *W,
                    DOUBLE           *WERR,
                    DOUBLE           *WGAP,
                    const INTEGER    *IBLOCK,
                    const INTEGER    *INDEXW,
                    const DOUBLE     *GERS,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *ISUPPZ,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- zlarscl2 ------------------------------------------------------------------
void
LAPACK_IMPL(zlarscl2)(const INTEGER        *M,
                      const INTEGER        *N,
                      const DOUBLE         *D,
                      DOUBLE_COMPLEX       *X,
                      const INTEGER        *LDX);

//-- zlartg --------------------------------------------------------------------
void
LAPACK_IMPL(zlartg)(const DOUBLE_COMPLEX     *F,
                    const DOUBLE_COMPLEX     *G,
                    DOUBLE                   *CS,
                    DOUBLE_COMPLEX           *SN,
                    DOUBLE_COMPLEX           *R);

//-- zlartv --------------------------------------------------------------------
void
LAPACK_IMPL(zlartv)(const INTEGER            *N,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *INCX,
                    DOUBLE_COMPLEX           *Y,
                    const INTEGER            *INCY,
                    const DOUBLE             *C,
                    const DOUBLE_COMPLEX     *S,
                    const INTEGER            *INCC);

//-- zlarz ---------------------------------------------------------------------
void
LAPACK_IMPL(zlarz)(const char               *SIDE,
                   const INTEGER            *M,
                   const INTEGER            *N,
                   const INTEGER            *L,
                   const DOUBLE_COMPLEX     *V,
                   const INTEGER            *INCV,
                   const DOUBLE_COMPLEX     *TAU,
                   DOUBLE_COMPLEX           *C,
                   const INTEGER            *LDC,
                   DOUBLE_COMPLEX           *WORK);

//-- zlarzb --------------------------------------------------------------------
void
LAPACK_IMPL(zlarzb)(const char               *SIDE,
                    const char               *TRANS,
                    const char               *DIRECT,
                    const char               *STOREV,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const INTEGER            *L,
                    const DOUBLE_COMPLEX     *V,
                    const INTEGER            *LDV,
                    const DOUBLE_COMPLEX     *T,
                    const INTEGER            *LDT,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LDWORK);

//-- zlarzt --------------------------------------------------------------------
void
LAPACK_IMPL(zlarzt)(const char               *DIRECT,
                    const char               *STOREV,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *V,
                    const INTEGER            *LDV,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *T,
                    const INTEGER            *LDT);

//-- zlascl --------------------------------------------------------------------
void
LAPACK_IMPL(zlascl)(const char       *TYPE,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    const DOUBLE     *CFROM,
                    const DOUBLE     *CTO,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- zlascl2 -------------------------------------------------------------------
void
LAPACK_IMPL(zlascl2)(const INTEGER    *M,
                     const INTEGER    *N,
                     const DOUBLE     *D,
                     DOUBLE_COMPLEX   *X,
                     const INTEGER    *LDX);

//-- zlaset --------------------------------------------------------------------
void
LAPACK_IMPL(zlaset)(const char               *UPLO,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *ALPHA,
                    const DOUBLE_COMPLEX     *BETA,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA);

//-- zlasr ---------------------------------------------------------------------
void
LAPACK_IMPL(zlasr)(const char           *SIDE,
                   const char           *PIVOT,
                   const char           *DIRECT,
                   const INTEGER        *M,
                   const INTEGER        *N,
                   const DOUBLE         *C,
                   const DOUBLE         *S,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA);

//-- zlassq --------------------------------------------------------------------
void
LAPACK_IMPL(zlassq)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *X,
                    const INTEGER            *INCX,
                    DOUBLE                   *SCALE,
                    DOUBLE                   *SUMSQ);

//-- zlaswp --------------------------------------------------------------------
void
LAPACK_IMPL(zlaswp)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const INTEGER    *K1,
                    const INTEGER    *K2,
                    const INTEGER    *IPIV,
                    const INTEGER    *INCX);

//-- zlasyf --------------------------------------------------------------------
void
LAPACK_IMPL(zlasyf)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NB,
                    INTEGER          *KB,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    DOUBLE_COMPLEX   *W,
                    const INTEGER    *LDW,
                    INTEGER          *INFO);

//-- zlat2c --------------------------------------------------------------------
void
LAPACK_IMPL(zlat2c)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    FLOAT_COMPLEX            *SA,
                    const INTEGER            *LDSA,
                    INTEGER                  *INFO);

//-- zlatbs --------------------------------------------------------------------
void
LAPACK_IMPL(zlatbs)(const char               *UPLO,
                    const char               *TRANS,
                    const char               *DIAG,
                    const char               *NORMIN,
                    const INTEGER            *N,
                    const INTEGER            *KD,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE_COMPLEX           *X,
                    DOUBLE                   *SCALE,
                    DOUBLE                   *CNORM,
                    INTEGER                  *INFO);

//-- zlatdf --------------------------------------------------------------------
void
LAPACK_IMPL(zlatdf)(const INTEGER            *IJOB,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *Z,
                    const INTEGER            *LDZ,
                    DOUBLE_COMPLEX           *RHS,
                    DOUBLE                   *RDSUM,
                    DOUBLE                   *RDSCAL,
                    const INTEGER            *IPIV,
                    const INTEGER            *JPIV);

//-- zlatps --------------------------------------------------------------------
void
LAPACK_IMPL(zlatps)(const char               *UPLO,
                    const char               *TRANS,
                    const char               *DIAG,
                    const char               *NORMIN,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE_COMPLEX           *X,
                    DOUBLE                   *SCALE,
                    DOUBLE                   *CNORM,
                    INTEGER                  *INFO);

//-- zlatrd --------------------------------------------------------------------
void
LAPACK_IMPL(zlatrd)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NB,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *W,
                    const INTEGER    *LDW);

//-- zlatrs --------------------------------------------------------------------
void
LAPACK_IMPL(zlatrs)(const char               *UPLO,
                    const char               *TRANS,
                    const char               *DIAG,
                    const char               *NORMIN,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *X,
                    DOUBLE                   *SCALE,
                    DOUBLE                   *CNORM,
                    INTEGER                  *INFO);

//-- zlatrz --------------------------------------------------------------------
void
LAPACK_IMPL(zlatrz)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *L,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK);

//-- zlatzm --------------------------------------------------------------------
void
LAPACK_IMPL(zlatzm)(const char               *SIDE,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *V,
                    const INTEGER            *INCV,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C1,
                    DOUBLE_COMPLEX           *C2,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK);

//-- zlauu2 --------------------------------------------------------------------
void
LAPACK_IMPL(zlauu2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- zlauum --------------------------------------------------------------------
void
LAPACK_IMPL(zlauum)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- zpbcon --------------------------------------------------------------------
void
LAPACK_IMPL(zpbcon)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *KD,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zpbequ --------------------------------------------------------------------
void
LAPACK_IMPL(zpbequ)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *KD,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE                   *S,
                    DOUBLE                   *SCOND,
                    DOUBLE                   *AMAX,
                    INTEGER                  *INFO);

//-- zpbrfs --------------------------------------------------------------------
void
LAPACK_IMPL(zpbrfs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *KD,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    const DOUBLE_COMPLEX     *AFB,
                    const INTEGER            *LDAFB,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zpbstf --------------------------------------------------------------------
void
LAPACK_IMPL(zpbstf)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO);

//-- zpbsv ---------------------------------------------------------------------
void
LAPACK_IMPL(zpbsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *KD,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *AB,
                   const INTEGER        *LDAB,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zpbsvx --------------------------------------------------------------------
void
LAPACK_IMPL(zpbsvx)(const char       *FACT,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    const INTEGER    *NRHS,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    DOUBLE_COMPLEX   *AFB,
                    const INTEGER    *LDAFB,
                    char             *EQUED,
                    DOUBLE           *S,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *LDX,
                    DOUBLE           *RCOND,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zpbtf2 --------------------------------------------------------------------
void
LAPACK_IMPL(zpbtf2)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO);

//-- zpbtrf --------------------------------------------------------------------
void
LAPACK_IMPL(zpbtrf)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO);

//-- zpbtrs --------------------------------------------------------------------
void
LAPACK_IMPL(zpbtrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *KD,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zpftrf --------------------------------------------------------------------
void
LAPACK_IMPL(zpftrf)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    INTEGER          *INFO);

//-- zpftri --------------------------------------------------------------------
void
LAPACK_IMPL(zpftri)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    INTEGER          *INFO);

//-- zpftrs --------------------------------------------------------------------
void
LAPACK_IMPL(zpftrs)(const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zpocon --------------------------------------------------------------------
void
LAPACK_IMPL(zpocon)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zpoequ --------------------------------------------------------------------
void
LAPACK_IMPL(zpoequ)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *S,
                    DOUBLE                   *SCOND,
                    DOUBLE                   *AMAX,
                    INTEGER                  *INFO);

//-- zpoequb -------------------------------------------------------------------
void
LAPACK_IMPL(zpoequb)(const INTEGER            *N,
                     const DOUBLE_COMPLEX     *A,
                     const INTEGER            *LDA,
                     DOUBLE                   *S,
                     DOUBLE                   *SCOND,
                     DOUBLE                   *AMAX,
                     INTEGER                  *INFO);

//-- zporfs --------------------------------------------------------------------
void
LAPACK_IMPL(zporfs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *AF,
                    const INTEGER            *LDAF,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zporfsx -------------------------------------------------------------------
/*
void
LAPACK_IMPL(zporfsx)(const char               *UPLO,
                     const char               *EQUED,
                     const INTEGER            *N,
                     const INTEGER            *NRHS,
                     const DOUBLE_COMPLEX     *A,
                     const INTEGER            *LDA,
                     const DOUBLE_COMPLEX     *AF,
                     const INTEGER            *LDAF,
                     DOUBLE                   *S,
                     const DOUBLE_COMPLEX     *B,
                     const INTEGER            *LDB,
                     DOUBLE_COMPLEX           *X,
                     const INTEGER            *LDX,
                     DOUBLE                   *RCOND,
                     DOUBLE                   *BERR,
                     const INTEGER            *N_ERR_BNDS,
                     DOUBLE                   *ERR_BNDS_NORM,
                     DOUBLE                   *ERR_BNDS_COMP,
                     const INTEGER            *NPARAMS,
                     DOUBLE                   *PARAMS,
                     DOUBLE_COMPLEX           *WORK,
                     DOUBLE                   *RWORK,
                     INTEGER                  *INFO);
*/

//-- zposv ---------------------------------------------------------------------
void
LAPACK_IMPL(zposv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zposvx --------------------------------------------------------------------
void
LAPACK_IMPL(zposvx)(const char       *FACT,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *AF,
                    const INTEGER    *LDAF,
                    char             *EQUED,
                    DOUBLE           *S,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *LDX,
                    DOUBLE           *RCOND,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zposvxx -------------------------------------------------------------------
/*
void
LAPACK_IMPL(zposvxx)(const char       *FACT,
                     const char       *UPLO,
                     const INTEGER    *N,
                     const INTEGER    *NRHS,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     DOUBLE_COMPLEX   *AF,
                     const INTEGER    *LDAF,
                     char             *EQUED,
                     DOUBLE           *S,
                     DOUBLE_COMPLEX   *B,
                     const INTEGER    *LDB,
                     DOUBLE_COMPLEX   *X,
                     const INTEGER    *LDX,
                     DOUBLE           *RCOND,
                     DOUBLE           *RPVGRW,
                     DOUBLE           *BERR,
                     const INTEGER    *N_ERR_BNDS,
                     DOUBLE           *ERR_BNDS_NORM,
                     DOUBLE           *ERR_BNDS_COMP,
                     const INTEGER    *NPARAMS,
                     DOUBLE           *PARAMS,
                     DOUBLE_COMPLEX   *WORK,
                     DOUBLE           *RWORK,
                     INTEGER          *INFO);
*/

//-- zpotf2 --------------------------------------------------------------------
void
LAPACK_IMPL(zpotf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- zpotrf --------------------------------------------------------------------
void
LAPACK_IMPL(zpotrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- zpotri --------------------------------------------------------------------
void
LAPACK_IMPL(zpotri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- zpotrs --------------------------------------------------------------------
void
LAPACK_IMPL(zpotrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zppcon --------------------------------------------------------------------
void
LAPACK_IMPL(zppcon)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zppequ --------------------------------------------------------------------
void
LAPACK_IMPL(zppequ)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE                   *S,
                    DOUBLE                   *SCOND,
                    DOUBLE                   *AMAX,
                    INTEGER                  *INFO);

//-- zpprfs --------------------------------------------------------------------
void
LAPACK_IMPL(zpprfs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    const DOUBLE_COMPLEX     *AFP,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zppsv ---------------------------------------------------------------------
void
LAPACK_IMPL(zppsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *AP,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zppsvx --------------------------------------------------------------------
void
LAPACK_IMPL(zppsvx)(const char       *FACT,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE_COMPLEX   *AP,
                    DOUBLE_COMPLEX   *AFP,
                    char             *EQUED,
                    DOUBLE           *S,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *LDX,
                    DOUBLE           *RCOND,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- zpptrf --------------------------------------------------------------------
void
LAPACK_IMPL(zpptrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    INTEGER          *INFO);

//-- zpptri --------------------------------------------------------------------
void
LAPACK_IMPL(zpptri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    INTEGER          *INFO);

//-- zpptrs --------------------------------------------------------------------
void
LAPACK_IMPL(zpptrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zpstf2 --------------------------------------------------------------------
void
LAPACK_IMPL(zpstf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *PIV,
                    INTEGER          *RANK,
                    const DOUBLE     *TOL,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- zpstrf --------------------------------------------------------------------
void
LAPACK_IMPL(zpstrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *PIV,
                    INTEGER          *RANK,
                    const DOUBLE     *TOL,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- zptcon --------------------------------------------------------------------
void
LAPACK_IMPL(zptcon)(const INTEGER            *N,
                    const DOUBLE             *D,
                    const DOUBLE_COMPLEX     *E,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zpteqr --------------------------------------------------------------------
void
LAPACK_IMPL(zpteqr)(const char       *COMPZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- zptrfs --------------------------------------------------------------------
void
LAPACK_IMPL(zptrfs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE             *D,
                    const DOUBLE_COMPLEX     *E,
                    const DOUBLE             *DF,
                    const DOUBLE_COMPLEX     *EF,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zptsv ---------------------------------------------------------------------
void
LAPACK_IMPL(zptsv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *D,
                   DOUBLE_COMPLEX       *E,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zptsvx --------------------------------------------------------------------
void
LAPACK_IMPL(zptsvx)(const char               *FACT,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE             *D,
                    const DOUBLE_COMPLEX     *E,
                    DOUBLE                   *DF,
                    DOUBLE_COMPLEX           *EF,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *RCOND,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zpttrf --------------------------------------------------------------------
void
LAPACK_IMPL(zpttrf)(const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE_COMPLEX   *E,
                    INTEGER          *INFO);

//-- zpttrs --------------------------------------------------------------------
void
LAPACK_IMPL(zpttrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE             *D,
                    const DOUBLE_COMPLEX     *E,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zptts2 --------------------------------------------------------------------
void
LAPACK_IMPL(zptts2)(const INTEGER            *IUPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE             *D,
                    const DOUBLE_COMPLEX     *E,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB);

//-- zrot ----------------------------------------------------------------------
void
LAPACK_IMPL(zrot)(const INTEGER            *N,
                  DOUBLE_COMPLEX           *CX,
                  const INTEGER            *INCX,
                  DOUBLE_COMPLEX           *CY,
                  const INTEGER            *INCY,
                  const DOUBLE             *C,
                  const DOUBLE_COMPLEX     *S);

//-- zspcon --------------------------------------------------------------------
void
LAPACK_IMPL(zspcon)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    const INTEGER            *IPIV,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zspmv ---------------------------------------------------------------------
void
LAPACK_IMPL(zspmv)(const char               *UPLO,
                   const INTEGER            *N,
                   const DOUBLE_COMPLEX     *ALPHA,
                   const DOUBLE_COMPLEX     *AP,
                   const DOUBLE_COMPLEX     *X,
                   const INTEGER            *INCX,
                   const DOUBLE_COMPLEX     *BETA,
                   DOUBLE_COMPLEX           *Y,
                   const INTEGER            *INCY);

//-- zspr ----------------------------------------------------------------------
void
LAPACK_IMPL(zspr)(const char               *UPLO,
                  const INTEGER            *N,
                  const DOUBLE_COMPLEX     *ALPHA,
                  const DOUBLE_COMPLEX     *X,
                  const INTEGER            *INCX,
                  DOUBLE_COMPLEX           *AP);

//-- zsprfs --------------------------------------------------------------------
void
LAPACK_IMPL(zsprfs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    const DOUBLE_COMPLEX     *AFP,
                    const INTEGER            *IPIV,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zspsv ---------------------------------------------------------------------
void
LAPACK_IMPL(zspsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *AP,
                   INTEGER              *IPIV,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zspsvx --------------------------------------------------------------------
void
LAPACK_IMPL(zspsvx)(const char               *FACT,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE_COMPLEX           *AFP,
                    INTEGER                  *IPIV,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *RCOND,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zsptrf --------------------------------------------------------------------
void
LAPACK_IMPL(zsptrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zsptri --------------------------------------------------------------------
void
LAPACK_IMPL(zsptri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    const INTEGER    *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zsptrs --------------------------------------------------------------------
void
LAPACK_IMPL(zsptrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    const INTEGER            *IPIV,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zstedc --------------------------------------------------------------------
void
LAPACK_IMPL(zstedc)(const char       *COMPZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    const INTEGER    *LRWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- zstegr --------------------------------------------------------------------
void
LAPACK_IMPL(zstegr)(const char       *JOBZ,
                    const char       *RANGE,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    const DOUBLE     *ABSTOL,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *ISUPPZ,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- zstein --------------------------------------------------------------------
void
LAPACK_IMPL(zstein)(const INTEGER    *N,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    const INTEGER    *M,
                    const DOUBLE     *W,
                    const INTEGER    *IBLOCK,
                    const INTEGER    *ISPLIT,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *IFAIL,
                    INTEGER          *INFO);

//-- zstemr --------------------------------------------------------------------
void
LAPACK_IMPL(zstemr)(const char       *JOBZ,
                    const char       *RANGE,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    const DOUBLE     *VL,
                    const DOUBLE     *VU,
                    const INTEGER    *IL,
                    const INTEGER    *IU,
                    INTEGER          *M,
                    DOUBLE           *W,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    const INTEGER    *NZC,
                    INTEGER          *ISUPPZ,
                    LOGICAL          *TRYRAC,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- zsteqr --------------------------------------------------------------------
void
LAPACK_IMPL(zsteqr)(const char       *COMPZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- zsycon --------------------------------------------------------------------
void
LAPACK_IMPL(zsycon)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const INTEGER            *IPIV,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zsyconv -------------------------------------------------------------------
void
LAPACK_IMPL(zsyconv)(const char       *UPLO,
                     const char       *WAY,
                     const INTEGER    *N,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE_COMPLEX   *WORK,
                     INTEGER          *INFO);

//-- zsyequb -------------------------------------------------------------------
void
LAPACK_IMPL(zsyequb)(const char               *UPLO,
                     const INTEGER            *N,
                     const DOUBLE_COMPLEX     *A,
                     const INTEGER            *LDA,
                     DOUBLE                   *S,
                     DOUBLE                   *SCOND,
                     DOUBLE                   *AMAX,
                     DOUBLE_COMPLEX           *WORK,
                     INTEGER                  *INFO);

//-- zsymv ---------------------------------------------------------------------
void
LAPACK_IMPL(zsymv)(const char               *UPLO,
                   const INTEGER            *N,
                   const DOUBLE_COMPLEX     *ALPHA,
                   const DOUBLE_COMPLEX     *A,
                   const INTEGER            *LDA,
                   const DOUBLE_COMPLEX     *X,
                   const INTEGER            *INCX,
                   const DOUBLE_COMPLEX     *BETA,
                   DOUBLE_COMPLEX           *Y,
                   const INTEGER            *INCY);

//-- zsyr ----------------------------------------------------------------------
void
LAPACK_IMPL(zsyr)(const char               *UPLO,
                  const INTEGER            *N,
                  const DOUBLE_COMPLEX     *ALPHA,
                  const DOUBLE_COMPLEX     *X,
                  const INTEGER            *INCX,
                  DOUBLE_COMPLEX           *A,
                  const INTEGER            *LDA);

//-- zsyrfs --------------------------------------------------------------------
void
LAPACK_IMPL(zsyrfs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *AF,
                    const INTEGER            *LDAF,
                    const INTEGER            *IPIV,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zsyrfsx -------------------------------------------------------------------
/*
void
LAPACK_IMPL(zsyrfsx)(const char               *UPLO,
                     const char               *EQUED,
                     const INTEGER            *N,
                     const INTEGER            *NRHS,
                     const DOUBLE_COMPLEX     *A,
                     const INTEGER            *LDA,
                     const DOUBLE_COMPLEX     *AF,
                     const INTEGER            *LDAF,
                     const INTEGER            *IPIV,
                     DOUBLE                   *S,
                     const DOUBLE_COMPLEX     *B,
                     const INTEGER            *LDB,
                     DOUBLE_COMPLEX           *X,
                     const INTEGER            *LDX,
                     DOUBLE                   *RCOND,
                     DOUBLE                   *BERR,
                     const INTEGER            *N_ERR_BNDS,
                     DOUBLE                   *ERR_BNDS_NORM,
                     DOUBLE                   *ERR_BNDS_COMP,
                     const INTEGER            *NPARAMS,
                     DOUBLE                   *PARAMS,
                     DOUBLE_COMPLEX           *WORK,
                     DOUBLE                   *RWORK,
                     INTEGER                  *INFO);
*/

//-- zsysv ---------------------------------------------------------------------
void
LAPACK_IMPL(zsysv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   INTEGER              *IPIV,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   DOUBLE_COMPLEX       *WORK,
                   const INTEGER        *LWORK,
                   INTEGER              *INFO);

//-- zsysvx --------------------------------------------------------------------
void
LAPACK_IMPL(zsysvx)(const char               *FACT,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *AF,
                    const INTEGER            *LDAF,
                    INTEGER                  *IPIV,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *RCOND,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zsysvxx -------------------------------------------------------------------
/*
void
LAPACK_IMPL(zsysvxx)(const char       *FACT,
                     const char       *UPLO,
                     const INTEGER    *N,
                     const INTEGER    *NRHS,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     DOUBLE_COMPLEX   *AF,
                     const INTEGER    *LDAF,
                     INTEGER          *IPIV,
                     char             *EQUED,
                     DOUBLE           *S,
                     DOUBLE_COMPLEX   *B,
                     const INTEGER    *LDB,
                     DOUBLE_COMPLEX   *X,
                     const INTEGER    *LDX,
                     DOUBLE           *RCOND,
                     DOUBLE           *RPVGRW,
                     DOUBLE           *BERR,
                     const INTEGER    *N_ERR_BNDS,
                     DOUBLE           *ERR_BNDS_NORM,
                     DOUBLE           *ERR_BNDS_COMP,
                     const INTEGER    *NPARAMS,
                     DOUBLE           *PARAMS,
                     DOUBLE_COMPLEX   *WORK,
                     DOUBLE           *RWORK,
                     INTEGER          *INFO);
*/

//-- zsyswapr ------------------------------------------------------------------
void
LAPACK_IMPL(zsyswapr)(const char       *UPLO,
                      const INTEGER    *N,
                      DOUBLE_COMPLEX   *A,
                      const INTEGER    *LDA,
                      const INTEGER    *I1,
                      const INTEGER    *I2);

//-- zsytf2 --------------------------------------------------------------------
void
LAPACK_IMPL(zsytf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zsytrf --------------------------------------------------------------------
void
LAPACK_IMPL(zsytrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zsytri --------------------------------------------------------------------
void
LAPACK_IMPL(zsytri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zsytri2 -------------------------------------------------------------------
void
LAPACK_IMPL(zsytri2)(const char       *UPLO,
                     const INTEGER    *N,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE_COMPLEX   *WORK,
                     const INTEGER    *LWORK,
                     INTEGER          *INFO);

//-- zsytri2x ------------------------------------------------------------------
void
LAPACK_IMPL(zsytri2x)(const char       *UPLO,
                      const INTEGER    *N,
                      DOUBLE_COMPLEX   *A,
                      const INTEGER    *LDA,
                      const INTEGER    *IPIV,
                      DOUBLE_COMPLEX   *WORK,
                      const INTEGER    *NB,
                      INTEGER          *INFO);

//-- zsytrs --------------------------------------------------------------------
void
LAPACK_IMPL(zsytrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const INTEGER            *IPIV,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zsytrs2 -------------------------------------------------------------------
void
LAPACK_IMPL(zsytrs2)(const char       *UPLO,
                     const INTEGER    *N,
                     const INTEGER    *NRHS,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE_COMPLEX   *B,
                     const INTEGER    *LDB,
                     DOUBLE_COMPLEX   *WORK,
                     INTEGER          *INFO);

//-- ztbcon --------------------------------------------------------------------
void
LAPACK_IMPL(ztbcon)(const char               *NORM,
                    const char               *UPLO,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const INTEGER            *KD,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- ztbrfs --------------------------------------------------------------------
void
LAPACK_IMPL(ztbrfs)(const char               *UPLO,
                    const char               *TRANS,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const INTEGER            *KD,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    const DOUBLE_COMPLEX     *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- ztbtrs --------------------------------------------------------------------
void
LAPACK_IMPL(ztbtrs)(const char               *UPLO,
                    const char               *TRANS,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const INTEGER            *KD,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- ztfsm ---------------------------------------------------------------------
void
LAPACK_IMPL(ztfsm)(const char               *TRANSR,
                   const char               *SIDE,
                   const char               *UPLO,
                   const char               *TRANS,
                   const char               *DIAG,
                   const INTEGER            *M,
                   const INTEGER            *N,
                   const DOUBLE_COMPLEX     *ALPHA,
                   const DOUBLE_COMPLEX     *A,
                   DOUBLE_COMPLEX           *B,
                   const INTEGER            *LDB);

//-- ztftri --------------------------------------------------------------------
void
LAPACK_IMPL(ztftri)(const char       *TRANSR,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    INTEGER          *INFO);

//-- ztfttp --------------------------------------------------------------------
void
LAPACK_IMPL(ztfttp)(const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *ARF,
                    DOUBLE_COMPLEX           *AP,
                    INTEGER                  *INFO);

//-- ztfttr --------------------------------------------------------------------
void
LAPACK_IMPL(ztfttr)(const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *ARF,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    INTEGER                  *INFO);

//-- ztgevc --------------------------------------------------------------------
void
LAPACK_IMPL(ztgevc)(const char               *SIDE,
                    const char               *HOWMNY,
                    const LOGICAL            *SELECT,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *S,
                    const INTEGER            *LDS,
                    const DOUBLE_COMPLEX     *P,
                    const INTEGER            *LDP,
                    DOUBLE_COMPLEX           *VL,
                    const INTEGER            *LDVL,
                    DOUBLE_COMPLEX           *VR,
                    const INTEGER            *LDVR,
                    const INTEGER            *MM,
                    INTEGER                  *M,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- ztgex2 --------------------------------------------------------------------
void
LAPACK_IMPL(ztgex2)(const LOGICAL    *WANTQ,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    const INTEGER    *J1,
                    INTEGER          *INFO);

//-- ztgexc --------------------------------------------------------------------
void
LAPACK_IMPL(ztgexc)(const LOGICAL    *WANTQ,
                    const LOGICAL    *WANTZ,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    const INTEGER    *IFST,
                    INTEGER          *ILST,
                    INTEGER          *INFO);

//-- ztgsen --------------------------------------------------------------------
void
LAPACK_IMPL(ztgsen)(const INTEGER    *IJOB,
                    const LOGICAL    *WANTQ,
                    const LOGICAL    *WANTZ,
                    const LOGICAL    *SELECT,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    DOUBLE_COMPLEX   *ALPHA,
                    DOUBLE_COMPLEX   *BETA,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    INTEGER          *M,
                    DOUBLE           *PL,
                    DOUBLE           *PR,
                    DOUBLE           *DIF,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    INTEGER          *INFO);

//-- ztgsja --------------------------------------------------------------------
void
LAPACK_IMPL(ztgsja)(const char       *JOBU,
                    const char       *JOBV,
                    const char       *JOBQ,
                    const INTEGER    *M,
                    const INTEGER    *P,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const INTEGER    *L,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB,
                    const DOUBLE     *TOLA,
                    const DOUBLE     *TOLB,
                    DOUBLE           *ALPHA,
                    DOUBLE           *BETA,
                    DOUBLE_COMPLEX   *U,
                    const INTEGER    *LDU,
                    DOUBLE_COMPLEX   *V,
                    const INTEGER    *LDV,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *NCYCLE,
                    INTEGER          *INFO);

//-- ztgsna --------------------------------------------------------------------
void
LAPACK_IMPL(ztgsna)(const char               *JOB,
                    const char               *HOWMNY,
                    const LOGICAL            *SELECT,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    const DOUBLE_COMPLEX     *VL,
                    const INTEGER            *LDVL,
                    const DOUBLE_COMPLEX     *VR,
                    const INTEGER            *LDVR,
                    DOUBLE                   *S,
                    DOUBLE                   *DIF,
                    const INTEGER            *MM,
                    INTEGER                  *M,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *IWORK,
                    INTEGER                  *INFO);

//-- ztgsy2 --------------------------------------------------------------------
void
LAPACK_IMPL(ztgsy2)(const char               *TRANS,
                    const INTEGER            *IJOB,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    const DOUBLE_COMPLEX     *D,
                    const INTEGER            *LDD,
                    const DOUBLE_COMPLEX     *E,
                    const INTEGER            *LDE,
                    DOUBLE_COMPLEX           *F,
                    const INTEGER            *LDF,
                    DOUBLE                   *SCALE,
                    DOUBLE                   *RDSUM,
                    DOUBLE                   *RDSCAL,
                    INTEGER                  *INFO);

//-- ztgsyl --------------------------------------------------------------------
void
LAPACK_IMPL(ztgsyl)(const char               *TRANS,
                    const INTEGER            *IJOB,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    const DOUBLE_COMPLEX     *D,
                    const INTEGER            *LDD,
                    const DOUBLE_COMPLEX     *E,
                    const INTEGER            *LDE,
                    DOUBLE_COMPLEX           *F,
                    const INTEGER            *LDF,
                    DOUBLE                   *SCALE,
                    DOUBLE                   *DIF,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *IWORK,
                    INTEGER                  *INFO);

//-- ztpcon --------------------------------------------------------------------
void
LAPACK_IMPL(ztpcon)(const char               *NORM,
                    const char               *UPLO,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- ztprfs --------------------------------------------------------------------
void
LAPACK_IMPL(ztprfs)(const char               *UPLO,
                    const char               *TRANS,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    const DOUBLE_COMPLEX     *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- ztptri --------------------------------------------------------------------
void
LAPACK_IMPL(ztptri)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    INTEGER          *INFO);

//-- ztptrs --------------------------------------------------------------------
void
LAPACK_IMPL(ztptrs)(const char               *UPLO,
                    const char               *TRANS,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- ztpttf --------------------------------------------------------------------
void
LAPACK_IMPL(ztpttf)(const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE_COMPLEX           *ARF,
                    INTEGER                  *INFO);

//-- ztpttr --------------------------------------------------------------------
void
LAPACK_IMPL(ztpttr)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    INTEGER                  *INFO);

//-- ztrcon --------------------------------------------------------------------
void
LAPACK_IMPL(ztrcon)(const char               *NORM,
                    const char               *UPLO,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- ztrevc --------------------------------------------------------------------
void
LAPACK_IMPL(ztrevc)(const char       *SIDE,
                    const char       *HOWMNY,
                    const LOGICAL    *SELECT,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *T,
                    const INTEGER    *LDT,
                    DOUBLE_COMPLEX   *VL,
                    const INTEGER    *LDVL,
                    DOUBLE_COMPLEX   *VR,
                    const INTEGER    *LDVR,
                    const INTEGER    *MM,
                    INTEGER          *M,
                    DOUBLE_COMPLEX   *WORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO);

//-- ztrexc --------------------------------------------------------------------
void
LAPACK_IMPL(ztrexc)(const char       *COMPQ,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *T,
                    const INTEGER    *LDT,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    const INTEGER    *IFST,
                    const INTEGER    *ILST,
                    INTEGER          *INFO);

//-- ztrrfs --------------------------------------------------------------------
void
LAPACK_IMPL(ztrrfs)(const char               *UPLO,
                    const char               *TRANS,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    const DOUBLE_COMPLEX     *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- ztrsen --------------------------------------------------------------------
void
LAPACK_IMPL(ztrsen)(const char       *JOB,
                    const char       *COMPQ,
                    const LOGICAL    *SELECT,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *T,
                    const INTEGER    *LDT,
                    DOUBLE_COMPLEX   *Q,
                    const INTEGER    *LDQ,
                    DOUBLE_COMPLEX   *W,
                    INTEGER          *M,
                    DOUBLE           *S,
                    DOUBLE           *SEP,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- ztrsna --------------------------------------------------------------------
void
LAPACK_IMPL(ztrsna)(const char               *JOB,
                    const char               *HOWMNY,
                    const LOGICAL            *SELECT,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *T,
                    const INTEGER            *LDT,
                    const DOUBLE_COMPLEX     *VL,
                    const INTEGER            *LDVL,
                    const DOUBLE_COMPLEX     *VR,
                    const INTEGER            *LDVR,
                    DOUBLE                   *S,
                    DOUBLE                   *SEP,
                    const INTEGER            *MM,
                    INTEGER                  *M,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LDWORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- ztrsyl --------------------------------------------------------------------
void
LAPACK_IMPL(ztrsyl)(const char               *TRANA,
                    const char               *TRANB,
                    const INTEGER            *ISGN,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE                   *SCALE,
                    INTEGER                  *INFO);

//-- ztrti2 --------------------------------------------------------------------
void
LAPACK_IMPL(ztrti2)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- ztrtri --------------------------------------------------------------------
void
LAPACK_IMPL(ztrtri)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- ztrtrs --------------------------------------------------------------------
void
LAPACK_IMPL(ztrtrs)(const char               *UPLO,
                    const char               *TRANS,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- ztrttf --------------------------------------------------------------------
void
LAPACK_IMPL(ztrttf)(const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *ARF,
                    INTEGER                  *INFO);

//-- ztrttp --------------------------------------------------------------------
void
LAPACK_IMPL(ztrttp)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *AP,
                    INTEGER                  *INFO);

//-- ztzrqf --------------------------------------------------------------------
void
LAPACK_IMPL(ztzrqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    INTEGER          *INFO);

//-- ztzrzf --------------------------------------------------------------------
void
LAPACK_IMPL(ztzrzf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zunbdb --------------------------------------------------------------------
void
LAPACK_IMPL(zunbdb)(const char       *TRANS,
                    const char       *SIGNS,
                    const INTEGER    *M,
                    const INTEGER    *P,
                    const INTEGER    *Q,
                    DOUBLE_COMPLEX   *X11,
                    const INTEGER    *LDX11,
                    DOUBLE_COMPLEX   *X12,
                    const INTEGER    *LDX12,
                    DOUBLE_COMPLEX   *X21,
                    const INTEGER    *LDX21,
                    DOUBLE_COMPLEX   *X22,
                    const INTEGER    *LDX22,
                    DOUBLE           *THETA,
                    DOUBLE           *PHI,
                    DOUBLE_COMPLEX   *TAUP1,
                    DOUBLE_COMPLEX   *TAUP2,
                    DOUBLE_COMPLEX   *TAUQ1,
                    DOUBLE_COMPLEX   *TAUQ2,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zuncsd --------------------------------------------------------------------
void
LAPACK_IMPL(zuncsd)(const char               *JOBU1,
                    const char               *JOBU2,
                    const char               *JOBV1T,
                    const char               *JOBV2T,
                    const char               *TRANS,
                    const char               *SIGNS,
                    const INTEGER            *M,
                    const INTEGER            *P,
                    const INTEGER            *Q,
                    const DOUBLE_COMPLEX     *X11,
                    const INTEGER            *LDX11,
                    const DOUBLE_COMPLEX     *X12,
                    const INTEGER            *LDX12,
                    const DOUBLE_COMPLEX     *X21,
                    const INTEGER            *LDX21,
                    const DOUBLE_COMPLEX     *X22,
                    const INTEGER            *LDX22,
                    DOUBLE                   *THETA,
                    DOUBLE_COMPLEX           *U1,
                    const INTEGER            *LDU1,
                    DOUBLE_COMPLEX           *U2,
                    const INTEGER            *LDU2,
                    DOUBLE_COMPLEX           *V1T,
                    const INTEGER            *LDV1T,
                    DOUBLE_COMPLEX           *V2T,
                    const INTEGER            *LDV2T,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    DOUBLE                   *RWORK,
                    const INTEGER            *LRWORK,
                    INTEGER                  *IWORK,
                    INTEGER                  *INFO);

//-- zung2l --------------------------------------------------------------------
void
LAPACK_IMPL(zung2l)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zung2r --------------------------------------------------------------------
void
LAPACK_IMPL(zung2r)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zungbr --------------------------------------------------------------------
void
LAPACK_IMPL(zungbr)(const char               *VECT,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zunghr --------------------------------------------------------------------
void
LAPACK_IMPL(zunghr)(const INTEGER            *N,
                    const INTEGER            *ILO,
                    const INTEGER            *IHI,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zungl2 --------------------------------------------------------------------
void
LAPACK_IMPL(zungl2)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zunglq --------------------------------------------------------------------
void
LAPACK_IMPL(zunglq)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zungql --------------------------------------------------------------------
void
LAPACK_IMPL(zungql)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zungqr --------------------------------------------------------------------
void
LAPACK_IMPL(zungqr)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zungr2 --------------------------------------------------------------------
void
LAPACK_IMPL(zungr2)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zungrq --------------------------------------------------------------------
void
LAPACK_IMPL(zungrq)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zungtr --------------------------------------------------------------------
void
LAPACK_IMPL(zungtr)(const char               *UPLO,
                    const INTEGER            *N,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zunm2l --------------------------------------------------------------------
void
LAPACK_IMPL(zunm2l)(const char               *SIDE,
                    const char               *TRANS,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zunm2r --------------------------------------------------------------------
void
LAPACK_IMPL(zunm2r)(const char               *SIDE,
                    const char               *TRANS,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zunmbr --------------------------------------------------------------------
void
LAPACK_IMPL(zunmbr)(const char               *VECT,
                    const char               *SIDE,
                    const char               *TRANS,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zunmhr --------------------------------------------------------------------
void
LAPACK_IMPL(zunmhr)(const char               *SIDE,
                    const char               *TRANS,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *ILO,
                    const INTEGER            *IHI,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zunml2 --------------------------------------------------------------------
void
LAPACK_IMPL(zunml2)(const char               *SIDE,
                    const char               *TRANS,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zunmlq --------------------------------------------------------------------
void
LAPACK_IMPL(zunmlq)(const char               *SIDE,
                    const char               *TRANS,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zunmql --------------------------------------------------------------------
void
LAPACK_IMPL(zunmql)(const char               *SIDE,
                    const char               *TRANS,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zunmqr --------------------------------------------------------------------
void
LAPACK_IMPL(zunmqr)(const char               *SIDE,
                    const char               *TRANS,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zunmr2 --------------------------------------------------------------------
void
LAPACK_IMPL(zunmr2)(const char               *SIDE,
                    const char               *TRANS,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zunmr3 --------------------------------------------------------------------
void
LAPACK_IMPL(zunmr3)(const char               *SIDE,
                    const char               *TRANS,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const INTEGER            *L,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zunmrq --------------------------------------------------------------------
void
LAPACK_IMPL(zunmrq)(const char               *SIDE,
                    const char               *TRANS,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zunmrz --------------------------------------------------------------------
void
LAPACK_IMPL(zunmrz)(const char               *SIDE,
                    const char               *TRANS,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const INTEGER            *L,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zunmtr --------------------------------------------------------------------
void
LAPACK_IMPL(zunmtr)(const char               *SIDE,
                    const char               *UPLO,
                    const char               *TRANS,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zupgtr --------------------------------------------------------------------
void
LAPACK_IMPL(zupgtr)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *Q,
                    const INTEGER            *LDQ,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zupmtr --------------------------------------------------------------------
void
LAPACK_IMPL(zupmtr)(const char               *SIDE,
                    const char               *UPLO,
                    const char               *TRANS,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- cgetrf --------------------------------------------------------------------
void
LAPACK_IMPL(cgetrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    FLOAT_COMPLEX    *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- cgetrs --------------------------------------------------------------------
void
LAPACK_IMPL(cgetrs)(const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    FLOAT_COMPLEX    *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    FLOAT_COMPLEX    *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- clag2z --------------------------------------------------------------------
void
LAPACK_IMPL(clag2z)(const INTEGER    *M,
                    const INTEGER    *N,
                    FLOAT_COMPLEX    *SA,
                    const INTEGER    *LDSA,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- cpotrf --------------------------------------------------------------------
void
LAPACK_IMPL(cpotrf)(const char       *UPLO,
                    const INTEGER    *N,
                    FLOAT_COMPLEX    *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- cpotrs --------------------------------------------------------------------
void
LAPACK_IMPL(cpotrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    FLOAT_COMPLEX    *A,
                    const INTEGER    *LDA,
                    FLOAT_COMPLEX    *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- ilaprec -------------------------------------------------------------------
INTEGER
LAPACK_IMPL(ilaprec)(const char   *PREC);

//-- chla_transtype ------------------------------------------------------------
char
LAPACK_IMPL(chla_transtype)(const INTEGER    *TRANS);

//-- claswp --------------------------------------------------------------------
void
LAPACK_IMPL(claswp)(const INTEGER    *N,
                    FLOAT_COMPLEX    *A,
                    const INTEGER    *LDA,
                    const INTEGER    *K1,
                    const INTEGER    *K2,
                    const INTEGER    *IPIV,
                    const INTEGER    *INCX);

//-- izmax1 --------------------------------------------------------------------
INTEGER
LAPACK_IMPL(izmax1)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *CX,
                    const INTEGER            *INCX);

//-- ilazlc --------------------------------------------------------------------
INTEGER
LAPACK_IMPL(ilazlc)(const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA);

//-- ilazlr --------------------------------------------------------------------
INTEGER
LAPACK_IMPL(ilazlr)(const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA);

//-- cpotf2 --------------------------------------------------------------------
void
LAPACK_IMPL(cpotf2)(const char       *UPLO,
                    const INTEGER    *N,
                    FLOAT_COMPLEX    *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- clacgv --------------------------------------------------------------------
void
LAPACK_IMPL(clacgv)(const INTEGER    *N,
                    FLOAT_COMPLEX    *X,
                    const INTEGER    *INCX);

