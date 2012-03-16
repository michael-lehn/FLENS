#ifndef INTEGER
#define INTEGER int
#endif

#ifndef FLOAT
#define FLOAT float
#endif

#ifndef DOUBLE
#define DOUBLE double
#endif

#ifndef FLOAT_COMPLEX
#define FLOAT_COMPLEX double
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
LAPACK_DECL(dbbcsd)(const char       *JOBU1,
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
LAPACK_DECL(dbdsdc)(const char       *UPLO,
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
LAPACK_DECL(dbdsqr)(const char       *UPLO,
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
LAPACK_DECL(ddisna)(const char       *JOB,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *D,
                    DOUBLE           *SEP,
                    INTEGER          *INFO);

//-- dgbbrd --------------------------------------------------------------------
void
LAPACK_DECL(dgbbrd)(const char       *VECT,
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
LAPACK_DECL(dgbcon)(const char       *NORM,
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
LAPACK_DECL(dgbequ)(const INTEGER    *M,
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
LAPACK_DECL(dgbequb)(const INTEGER    *M,
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
LAPACK_DECL(dgbrfs)(const char       *TRANS,
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
LAPACK_DECL(dgbsv)(const INTEGER        *N,
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
LAPACK_DECL(dgbsvx)(const char       *FACT,
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
LAPACK_DECL(dgbtf2)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- dgbtrf --------------------------------------------------------------------
void
LAPACK_DECL(dgbtrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- dgbtrs --------------------------------------------------------------------
void
LAPACK_DECL(dgbtrs)(const char       *TRANS,
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
LAPACK_DECL(dgebak)(const char       *JOB,
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
LAPACK_DECL(dgebal)(const char       *JOB,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *ILO,
                    INTEGER          *IHI,
                    DOUBLE           *SCALE,
                    INTEGER          *INFO);

//-- dgebd2 --------------------------------------------------------------------
void
LAPACK_DECL(dgebd2)(const INTEGER    *M,
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
LAPACK_DECL(dgebrd)(const INTEGER    *M,
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
LAPACK_DECL(dgecon)(const char       *NORM,
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
LAPACK_DECL(dgeequ)(const INTEGER    *M,
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
LAPACK_DECL(dgeequb)(const INTEGER    *M,
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
LAPACK_DECL(dgees)(const char           *JOBVS,
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
LAPACK_DECL(dgeesx)(const char       *JOBVS,
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
LAPACK_DECL(dgeev)(const char           *JOBVL,
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
LAPACK_DECL(dgeevx)(const char       *BALANC,
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
LAPACK_DECL(dgegs)(const char           *JOBVSL,
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
LAPACK_DECL(dgegv)(const char           *JOBVL,
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
LAPACK_DECL(dgehd2)(const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgehrd --------------------------------------------------------------------
void
LAPACK_DECL(dgehrd)(const INTEGER    *N,
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
LAPACK_DECL(dgejsv)(const char       *JOBA,
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
LAPACK_DECL(dgelq2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgelqf --------------------------------------------------------------------
void
LAPACK_DECL(dgelqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgels ---------------------------------------------------------------------
void
LAPACK_DECL(dgels)(const char           *TRANS,
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
LAPACK_DECL(dgelsd)(const INTEGER    *M,
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
LAPACK_DECL(dgelss)(const INTEGER    *M,
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
LAPACK_DECL(dgelsx)(const INTEGER    *M,
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
LAPACK_DECL(dgelsy)(const INTEGER    *M,
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
LAPACK_DECL(dgeql2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgeqlf --------------------------------------------------------------------
void
LAPACK_DECL(dgeqlf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgeqp3 --------------------------------------------------------------------
void
LAPACK_DECL(dgeqp3)(const INTEGER    *M,
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
LAPACK_DECL(dgeqpf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *JPVT,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgeqr2 --------------------------------------------------------------------
void
LAPACK_DECL(dgeqr2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgeqr2p -------------------------------------------------------------------
void
LAPACK_DECL(dgeqr2p)(const INTEGER    *M,
                     const INTEGER    *N,
                     DOUBLE           *A,
                     const INTEGER    *LDA,
                     DOUBLE           *TAU,
                     DOUBLE           *WORK,
                     INTEGER          *INFO);

//-- dgeqrf --------------------------------------------------------------------
void
LAPACK_DECL(dgeqrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgeqrfp -------------------------------------------------------------------
void
LAPACK_DECL(dgeqrfp)(const INTEGER    *M,
                     const INTEGER    *N,
                     DOUBLE           *A,
                     const INTEGER    *LDA,
                     DOUBLE           *TAU,
                     DOUBLE           *WORK,
                     const INTEGER    *LWORK,
                     INTEGER          *INFO);

//-- dgerfs --------------------------------------------------------------------
void
LAPACK_DECL(dgerfs)(const char       *TRANS,
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
LAPACK_DECL(dgerq2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dgerqf --------------------------------------------------------------------
void
LAPACK_DECL(dgerqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgesc2 --------------------------------------------------------------------
void
LAPACK_DECL(dgesc2)(const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *RHS,
                    const INTEGER    *IPIV,
                    const INTEGER    *JPIV,
                    DOUBLE           *SCALE);

//-- dgesdd --------------------------------------------------------------------
void
LAPACK_DECL(dgesdd)(const char       *JOBZ,
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
LAPACK_DECL(dgesv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   INTEGER              *IPIV,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- dgesvd --------------------------------------------------------------------
void
LAPACK_DECL(dgesvd)(const char       *JOBU,
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
LAPACK_DECL(dgesvj)(const char       *JOBA,
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
LAPACK_DECL(dgesvx)(const char       *FACT,
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
LAPACK_DECL(dgetc2)(const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *JPIV,
                    INTEGER          *INFO);

//-- dgetf2 --------------------------------------------------------------------
void
LAPACK_DECL(dgetf2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- dgetrf --------------------------------------------------------------------
void
LAPACK_DECL(dgetrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- dgetri --------------------------------------------------------------------
void
LAPACK_DECL(dgetri)(const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dgetrs --------------------------------------------------------------------
void
LAPACK_DECL(dgetrs)(const char       *TRANS,
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
LAPACK_DECL(dggbak)(const char       *JOB,
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
LAPACK_DECL(dggbal)(const char       *JOB,
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
LAPACK_DECL(dgges)(const char           *JOBVSL,
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
LAPACK_DECL(dggesx)(const char       *JOBVSL,
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
LAPACK_DECL(dggev)(const char           *JOBVL,
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
LAPACK_DECL(dggevx)(const char       *BALANC,
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
LAPACK_DECL(dggglm)(const INTEGER    *N,
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
LAPACK_DECL(dgghrd)(const char       *COMPQ,
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
LAPACK_DECL(dgglse)(const INTEGER    *M,
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
LAPACK_DECL(dggqrf)(const INTEGER    *N,
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
LAPACK_DECL(dggrqf)(const INTEGER    *M,
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
LAPACK_DECL(dggsvd)(const char       *JOBU,
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
LAPACK_DECL(dggsvp)(const char       *JOBU,
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
LAPACK_DECL(dgsvj0)(const char       *JOBV,
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
LAPACK_DECL(dgsvj1)(const char       *JOBV,
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
LAPACK_DECL(dgtcon)(const char       *NORM,
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
LAPACK_DECL(dgtrfs)(const char       *TRANS,
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
LAPACK_DECL(dgtsv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *DL,
                   DOUBLE               *D,
                   DOUBLE               *DU,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- dgtsvx --------------------------------------------------------------------
void
LAPACK_DECL(dgtsvx)(const char       *FACT,
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
LAPACK_DECL(dgttrf)(const INTEGER    *N,
                    DOUBLE           *DL,
                    DOUBLE           *D,
                    DOUBLE           *DU,
                    DOUBLE           *DU2,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- dgttrs --------------------------------------------------------------------
void
LAPACK_DECL(dgttrs)(const char       *TRANS,
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
LAPACK_DECL(dgtts2)(const INTEGER    *ITRANS,
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
LAPACK_DECL(dhgeqz)(const char       *JOB,
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
LAPACK_DECL(dhsein)(const char       *SIDE,
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
LAPACK_DECL(dhseqr)(const char       *JOB,
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
LAPACK_DECL(disnan)(const DOUBLE     *DIN);

//-- dla_gbamv -----------------------------------------------------------------
void
LAPACK_DECL(dla_gbamv)(const INTEGER        *TRANS,
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
LAPACK_DECL(dla_gbrcond)(const char       *TRANS,
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
LAPACK_DECL(dla_gbrpvgrw)(const INTEGER    *N,
                          const INTEGER    *KL,
                          const INTEGER    *KU,
                          const INTEGER    *NCOLS,
                          const DOUBLE     *AB,
                          const INTEGER    *LDAB,
                          const DOUBLE     *AFB,
                          const INTEGER    *LDAFB);

//-- dla_geamv -----------------------------------------------------------------
void
LAPACK_DECL(dla_geamv)(const INTEGER        *TRANS,
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
LAPACK_DECL(dla_gercond)(const char       *TRANS,
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
LAPACK_DECL(dla_lin_berr)(const INTEGER    *N,
                          const INTEGER    *NZ,
                          const INTEGER    *NRHS,
                          const DOUBLE     *RES,
                          const DOUBLE     *AYB,
                          DOUBLE           *BERR);

//-- dla_porcond ---------------------------------------------------------------
DOUBLE
LAPACK_DECL(dla_porcond)(const char       *UPLO,
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
LAPACK_DECL(dla_porpvgrw)(const char       *UPLO,
                          const INTEGER    *NCOLS,
                          const DOUBLE     *A,
                          const INTEGER    *LDA,
                          const DOUBLE     *AF,
                          const INTEGER    *LDAF,
                          const DOUBLE     *WORK);

//-- dla_rpvgrw ----------------------------------------------------------------
DOUBLE
LAPACK_DECL(dla_rpvgrw)(const INTEGER    *N,
                        const INTEGER    *NCOLS,
                        const DOUBLE     *A,
                        const INTEGER    *LDA,
                        const DOUBLE     *AF,
                        const INTEGER    *LDAF);

//-- dla_syamv -----------------------------------------------------------------
void
LAPACK_DECL(dla_syamv)(const INTEGER        *UPLO,
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
LAPACK_DECL(dla_syrcond)(const char       *UPLO,
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
LAPACK_DECL(dla_syrpvgrw)(const char       *UPLO,
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
LAPACK_DECL(dla_wwaddw)(const INTEGER    *N,
                        DOUBLE           *X,
                        DOUBLE           *Y,
                        const DOUBLE     *W);

//-- dlabad --------------------------------------------------------------------
void
LAPACK_DECL(dlabad)(DOUBLE   *SMALL,
                    DOUBLE   *LARGE);

//-- dlabrd --------------------------------------------------------------------
void
LAPACK_DECL(dlabrd)(const INTEGER    *M,
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
LAPACK_DECL(dlacn2)(const INTEGER    *N,
                    DOUBLE           *V,
                    DOUBLE           *X,
                    INTEGER          *ISGN,
                    DOUBLE           *EST,
                    INTEGER          *KASE,
                    INTEGER          *ISAVE);

//-- dlacon --------------------------------------------------------------------
void
LAPACK_DECL(dlacon)(const INTEGER    *N,
                    DOUBLE           *V,
                    DOUBLE           *X,
                    INTEGER          *ISGN,
                    DOUBLE           *EST,
                    INTEGER          *KASE);

//-- dlacpy --------------------------------------------------------------------
void
LAPACK_DECL(dlacpy)(const char       *UPLO,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB);

//-- dladiv --------------------------------------------------------------------
void
LAPACK_DECL(dladiv)(const DOUBLE     *A,
                    const DOUBLE     *B,
                    const DOUBLE     *C,
                    const DOUBLE     *D,
                    DOUBLE           *P,
                    DOUBLE           *Q);

//-- dlae2 ---------------------------------------------------------------------
void
LAPACK_DECL(dlae2)(const DOUBLE     *A,
                   const DOUBLE     *B,
                   const DOUBLE     *C,
                   DOUBLE           *RT1,
                   DOUBLE           *RT2);

//-- dlaebz --------------------------------------------------------------------
void
LAPACK_DECL(dlaebz)(const INTEGER    *IJOB,
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
LAPACK_DECL(dlaed0)(const INTEGER    *ICOMPQ,
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
LAPACK_DECL(dlaed1)(const INTEGER    *N,
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
LAPACK_DECL(dlaed2)(INTEGER          *K,
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
LAPACK_DECL(dlaed3)(const INTEGER    *K,
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
LAPACK_DECL(dlaed4)(const INTEGER    *N,
                    const INTEGER    *I,
                    const DOUBLE     *D,
                    const DOUBLE     *Z,
                    DOUBLE           *DELTA,
                    const DOUBLE     *RHO,
                    DOUBLE           *DLAM,
                    INTEGER          *INFO);

//-- dlaed5 --------------------------------------------------------------------
void
LAPACK_DECL(dlaed5)(const INTEGER    *I,
                    const DOUBLE     *D,
                    const DOUBLE     *Z,
                    DOUBLE           *DELTA,
                    const DOUBLE     *RHO,
                    DOUBLE           *DLAM);

//-- dlaed6 --------------------------------------------------------------------
void
LAPACK_DECL(dlaed6)(const INTEGER    *KNITER,
                    const LOGICAL    *ORGATI,
                    const DOUBLE     *RHO,
                    const DOUBLE     *D,
                    const DOUBLE     *Z,
                    const DOUBLE     *FINIT,
                    DOUBLE           *TAU,
                    INTEGER          *INFO);

//-- dlaed7 --------------------------------------------------------------------
void
LAPACK_DECL(dlaed7)(const INTEGER    *ICOMPQ,
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
LAPACK_DECL(dlaed8)(const INTEGER    *ICOMPQ,
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
LAPACK_DECL(dlaed9)(const INTEGER    *K,
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
LAPACK_DECL(dlaeda)(const INTEGER    *N,
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
LAPACK_DECL(dlaein)(const LOGICAL    *RIGHTV,
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
LAPACK_DECL(dlaev2)(const DOUBLE     *A,
                    const DOUBLE     *B,
                    const DOUBLE     *C,
                    DOUBLE           *RT1,
                    DOUBLE           *RT2,
                    DOUBLE           *CS1,
                    DOUBLE           *SN1);

//-- dlaexc --------------------------------------------------------------------
void
LAPACK_DECL(dlaexc)(const LOGICAL    *WANTQ,
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
LAPACK_DECL(dlag2)(const DOUBLE         *A,
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
LAPACK_DECL(dlag2s)(const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    FLOAT            *SA,
                    const INTEGER    *LDSA,
                    INTEGER          *INFO);

//-- dlags2 --------------------------------------------------------------------
void
LAPACK_DECL(dlags2)(const LOGICAL    *UPPER,
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
LAPACK_DECL(dlagtf)(const INTEGER    *N,
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
LAPACK_DECL(dlagtm)(const char       *TRANS,
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
LAPACK_DECL(dlagts)(const INTEGER    *JOB,
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
LAPACK_DECL(dlagv2)(DOUBLE           *A,
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
LAPACK_DECL(dlahqr)(const LOGICAL    *WANTT,
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
LAPACK_DECL(dlahr2)(const INTEGER    *N,
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
LAPACK_DECL(dlahrd)(const INTEGER    *N,
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
LAPACK_DECL(dlaic1)(const INTEGER    *JOB,
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
LAPACK_DECL(dlaisnan)(const DOUBLE     *DIN1,
                      const DOUBLE     *DIN2);

//-- dlaln2 --------------------------------------------------------------------
void
LAPACK_DECL(dlaln2)(const LOGICAL    *LTRANS,
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
LAPACK_DECL(dlals0)(const INTEGER    *ICOMPQ,
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
LAPACK_DECL(dlalsa)(const INTEGER    *ICOMPQ,
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
LAPACK_DECL(dlalsd)(const char       *UPLO,
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
LAPACK_DECL(dlamch)(const char   *CMACH);

//-- dlamrg --------------------------------------------------------------------
void
LAPACK_DECL(dlamrg)(const INTEGER    *N1,
                    const INTEGER    *N2,
                    const DOUBLE     *A,
                    const INTEGER    *DTRD1,
                    const INTEGER    *DTRD2,
                    INTEGER          *INDEX);

//-- dlaneg --------------------------------------------------------------------
INTEGER
LAPACK_DECL(dlaneg)(const INTEGER    *N,
                    const DOUBLE     *D,
                    const DOUBLE     *LLD,
                    const DOUBLE     *SIGMA,
                    const DOUBLE     *PIVMIN,
                    const INTEGER    *R);

//-- dlangb --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlangb)(const char       *NORM,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *WORK);

//-- dlange --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlange)(const char       *NORM,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *WORK);

//-- dlangt --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlangt)(const char       *NORM,
                    const INTEGER    *N,
                    const DOUBLE     *DL,
                    const DOUBLE     *D,
                    const DOUBLE     *DU);

//-- dlanhs --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlanhs)(const char       *NORM,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *WORK);

//-- dlansb --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlansb)(const char       *NORM,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *WORK);

//-- dlansf --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlansf)(const char       *NORM,
                    const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    DOUBLE           *WORK);

//-- dlansp --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlansp)(const char       *NORM,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *WORK);

//-- dlanst --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlanst)(const char       *NORM,
                    const INTEGER    *N,
                    const DOUBLE     *D,
                    const DOUBLE     *E);

//-- dlansy --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlansy)(const char       *NORM,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *WORK);

//-- dlantb --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlantb)(const char       *NORM,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *WORK);

//-- dlantp --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlantp)(const char       *NORM,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *WORK);

//-- dlantr --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlantr)(const char       *NORM,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *WORK);

//-- dlanv2 --------------------------------------------------------------------
void
LAPACK_DECL(dlanv2)(DOUBLE   *A,
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
LAPACK_DECL(dlapll)(const INTEGER    *N,
                    DOUBLE           *X,
                    const INTEGER    *INCX,
                    DOUBLE           *Y,
                    const INTEGER    *INCY,
                    DOUBLE           *SSMIN);

//-- dlapmr --------------------------------------------------------------------
void
LAPACK_DECL(dlapmr)(const LOGICAL    *FORWRD,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    INTEGER          *K);

//-- dlapmt --------------------------------------------------------------------
void
LAPACK_DECL(dlapmt)(const LOGICAL    *FORWRD,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    INTEGER          *K);

//-- dlapy2 --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlapy2)(const DOUBLE     *X,
                    const DOUBLE     *Y);

//-- dlapy3 --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlapy3)(const DOUBLE     *X,
                    const DOUBLE     *Y,
                    const DOUBLE     *Z);

//-- dlaqgb --------------------------------------------------------------------
void
LAPACK_DECL(dlaqgb)(const INTEGER    *M,
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
LAPACK_DECL(dlaqge)(const INTEGER    *M,
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
LAPACK_DECL(dlaqp2)(const INTEGER    *M,
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
LAPACK_DECL(dlaqps)(const INTEGER    *M,
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
LAPACK_DECL(dlaqr0)(const LOGICAL    *WANTT,
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
LAPACK_DECL(dlaqr1)(const INTEGER    *N,
                    const DOUBLE     *H,
                    const INTEGER    *LDH,
                    const DOUBLE     *SR1,
                    const DOUBLE     *SI1,
                    const DOUBLE     *SR2,
                    const DOUBLE     *SI2,
                    DOUBLE           *V);

//-- dlaqr2 --------------------------------------------------------------------
void
LAPACK_DECL(dlaqr2)(const LOGICAL    *WANTT,
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
LAPACK_DECL(dlaqr3)(const LOGICAL    *WANTT,
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
LAPACK_DECL(dlaqr4)(const LOGICAL    *WANTT,
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
LAPACK_DECL(dlaqr5)(const LOGICAL    *WANTT,
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
LAPACK_DECL(dlaqsb)(const char       *UPLO,
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
LAPACK_DECL(dlaqsp)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- dlaqsy --------------------------------------------------------------------
void
LAPACK_DECL(dlaqsy)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- dlaqtr --------------------------------------------------------------------
void
LAPACK_DECL(dlaqtr)(const LOGICAL    *LTRAN,
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
LAPACK_DECL(dlar1v)(const INTEGER    *N,
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
LAPACK_DECL(dlar2v)(const INTEGER    *N,
                    DOUBLE           *X,
                    DOUBLE           *Y,
                    DOUBLE           *Z,
                    const INTEGER    *INCX,
                    const DOUBLE     *C,
                    const DOUBLE     *S,
                    const INTEGER    *INCC);

//-- dlarf ---------------------------------------------------------------------
void
LAPACK_DECL(dlarf)(const char           *SIDE,
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
LAPACK_DECL(dlarfb)(const char       *SIDE,
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
LAPACK_DECL(dlarfg)(const INTEGER    *N,
                    DOUBLE           *ALPHA,
                    DOUBLE           *X,
                    const INTEGER    *INCX,
                    DOUBLE           *TAU);

//-- dlarfgp -------------------------------------------------------------------
void
LAPACK_DECL(dlarfgp)(const INTEGER    *N,
                     DOUBLE           *ALPHA,
                     DOUBLE           *X,
                     const INTEGER    *INCX,
                     DOUBLE           *TAU);

//-- dlarft --------------------------------------------------------------------
void
LAPACK_DECL(dlarft)(const char       *DIRECT,
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
LAPACK_DECL(dlarfx)(const char       *SIDE,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *V,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK);

//-- dlargv --------------------------------------------------------------------
void
LAPACK_DECL(dlargv)(const INTEGER    *N,
                    DOUBLE           *X,
                    const INTEGER    *INCX,
                    DOUBLE           *Y,
                    const INTEGER    *INCY,
                    DOUBLE           *C,
                    const INTEGER    *INCC);

//-- dlarnv --------------------------------------------------------------------
void
LAPACK_DECL(dlarnv)(const INTEGER    *IDIST,
                    INTEGER          *ISEED,
                    const INTEGER    *N,
                    DOUBLE           *X);

//-- dlarra --------------------------------------------------------------------
void
LAPACK_DECL(dlarra)(const INTEGER    *N,
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
LAPACK_DECL(dlarrb)(const INTEGER    *N,
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
LAPACK_DECL(dlarrc)(const char       *JOBT,
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
LAPACK_DECL(dlarrd)(const char       *RANGE,
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
LAPACK_DECL(dlarre)(const char       *RANGE,
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
LAPACK_DECL(dlarrf)(const INTEGER    *N,
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
LAPACK_DECL(dlarrj)(const INTEGER    *N,
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
LAPACK_DECL(dlarrk)(const INTEGER    *N,
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
LAPACK_DECL(dlarrr)(const INTEGER    *N,
                    const DOUBLE     *D,
                    DOUBLE           *E,
                    INTEGER          *INFO);

//-- dlarrv --------------------------------------------------------------------
void
LAPACK_DECL(dlarrv)(const INTEGER    *N,
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
LAPACK_DECL(dlarscl2)(const INTEGER    *M,
                      const INTEGER    *N,
                      const DOUBLE     *D,
                      DOUBLE           *X,
                      const INTEGER    *LDX);

//-- dlartg --------------------------------------------------------------------
void
LAPACK_DECL(dlartg)(const DOUBLE     *F,
                    const DOUBLE     *G,
                    DOUBLE           *CS,
                    DOUBLE           *SN,
                    DOUBLE           *R);

//-- dlartgp -------------------------------------------------------------------
void
LAPACK_DECL(dlartgp)(const DOUBLE     *F,
                     const DOUBLE     *G,
                     DOUBLE           *CS,
                     DOUBLE           *SN,
                     DOUBLE           *R);

//-- dlartgs -------------------------------------------------------------------
void
LAPACK_DECL(dlartgs)(const DOUBLE     *X,
                     const DOUBLE     *Y,
                     const DOUBLE     *SIGMA,
                     DOUBLE           *CS,
                     DOUBLE           *SN);

//-- dlartv --------------------------------------------------------------------
void
LAPACK_DECL(dlartv)(const INTEGER    *N,
                    DOUBLE           *X,
                    const INTEGER    *INCX,
                    DOUBLE           *Y,
                    const INTEGER    *INCY,
                    const DOUBLE     *C,
                    const DOUBLE     *S,
                    const INTEGER    *INCC);

//-- dlaruv --------------------------------------------------------------------
void
LAPACK_DECL(dlaruv)(INTEGER          *ISEED,
                    const INTEGER    *N,
                    DOUBLE           *X);

//-- dlarz ---------------------------------------------------------------------
void
LAPACK_DECL(dlarz)(const char           *SIDE,
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
LAPACK_DECL(dlarzb)(const char       *SIDE,
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
LAPACK_DECL(dlarzt)(const char       *DIRECT,
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
LAPACK_DECL(dlas2)(const DOUBLE     *F,
                   const DOUBLE     *G,
                   const DOUBLE     *H,
                   DOUBLE           *SSMIN,
                   DOUBLE           *SSMAX);

//-- dlascl --------------------------------------------------------------------
void
LAPACK_DECL(dlascl)(const char       *TYPE,
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
LAPACK_DECL(dlascl2)(const INTEGER    *M,
                     const INTEGER    *N,
                     const DOUBLE     *D,
                     DOUBLE           *X,
                     const INTEGER    *LDX);

//-- dlasd0 --------------------------------------------------------------------
void
LAPACK_DECL(dlasd0)(const INTEGER    *N,
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
LAPACK_DECL(dlasd1)(const INTEGER    *NL,
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
LAPACK_DECL(dlasd2)(const INTEGER    *NL,
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
LAPACK_DECL(dlasd3)(const INTEGER    *NL,
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
LAPACK_DECL(dlasd4)(const INTEGER    *N,
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
LAPACK_DECL(dlasd5)(const INTEGER    *I,
                    const DOUBLE     *D,
                    const DOUBLE     *Z,
                    DOUBLE           *DELTA,
                    const DOUBLE     *RHO,
                    DOUBLE           *DSIGMA,
                    DOUBLE           *WORK);

//-- dlasd6 --------------------------------------------------------------------
void
LAPACK_DECL(dlasd6)(const INTEGER    *ICOMPQ,
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
LAPACK_DECL(dlasd7)(const INTEGER    *ICOMPQ,
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
LAPACK_DECL(dlasd8)(const INTEGER    *ICOMPQ,
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
LAPACK_DECL(dlasda)(const INTEGER    *ICOMPQ,
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
LAPACK_DECL(dlasdq)(const char       *UPLO,
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
LAPACK_DECL(dlasdt)(const INTEGER    *N,
                    INTEGER          *LVL,
                    INTEGER          *ND,
                    INTEGER          *INODE,
                    INTEGER          *NDIML,
                    INTEGER          *NDIMR,
                    const INTEGER    *MSUB);

//-- dlaset --------------------------------------------------------------------
void
LAPACK_DECL(dlaset)(const char       *UPLO,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *ALPHA,
                    const DOUBLE     *BETA,
                    DOUBLE           *A,
                    const INTEGER    *LDA);

//-- dlasq1 --------------------------------------------------------------------
void
LAPACK_DECL(dlasq1)(const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dlasq2 --------------------------------------------------------------------
void
LAPACK_DECL(dlasq2)(const INTEGER    *N,
                    DOUBLE           *Z,
                    INTEGER          *INFO);

//-- dlasq3 --------------------------------------------------------------------
void
LAPACK_DECL(dlasq3)(const INTEGER    *I0,
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
LAPACK_DECL(dlasq4)(const INTEGER    *I0,
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
LAPACK_DECL(dlasq5)(const INTEGER    *I0,
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
LAPACK_DECL(dlasq6)(const INTEGER    *I0,
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
LAPACK_DECL(dlasr)(const char           *SIDE,
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
LAPACK_DECL(dlasrt)(const char       *ID,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    INTEGER          *INFO);

//-- dlassq --------------------------------------------------------------------
void
LAPACK_DECL(dlassq)(const INTEGER    *N,
                    const DOUBLE     *X,
                    const INTEGER    *INCX,
                    DOUBLE           *SCALE,
                    DOUBLE           *SUMSQ);

//-- dlasv2 --------------------------------------------------------------------
void
LAPACK_DECL(dlasv2)(const DOUBLE     *F,
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
LAPACK_DECL(dlaswp)(const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const INTEGER    *K1,
                    const INTEGER    *K2,
                    const INTEGER    *IPIV,
                    const INTEGER    *INCX);

//-- dlasy2 --------------------------------------------------------------------
void
LAPACK_DECL(dlasy2)(const LOGICAL    *LTRANL,
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
LAPACK_DECL(dlasyf)(const char       *UPLO,
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
LAPACK_DECL(dlat2s)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    FLOAT            *SA,
                    const INTEGER    *LDSA,
                    INTEGER          *INFO);

//-- dlatbs --------------------------------------------------------------------
void
LAPACK_DECL(dlatbs)(const char       *UPLO,
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
LAPACK_DECL(dlatdf)(const INTEGER    *IJOB,
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
LAPACK_DECL(dlatps)(const char       *UPLO,
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
LAPACK_DECL(dlatrd)(const char       *UPLO,
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
LAPACK_DECL(dlatrs)(const char       *UPLO,
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
LAPACK_DECL(dlatrz)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *L,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK);

//-- dlatzm --------------------------------------------------------------------
void
LAPACK_DECL(dlatzm)(const char       *SIDE,
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
LAPACK_DECL(dlauu2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dlauum --------------------------------------------------------------------
void
LAPACK_DECL(dlauum)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dopgtr --------------------------------------------------------------------
void
LAPACK_DECL(dopgtr)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    const DOUBLE     *TAU,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dopmtr --------------------------------------------------------------------
void
LAPACK_DECL(dopmtr)(const char       *SIDE,
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
LAPACK_DECL(dorbdb)(const char       *TRANS,
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
LAPACK_DECL(dorcsd)(const char       *JOBU1,
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
LAPACK_DECL(dorg2l)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dorg2r --------------------------------------------------------------------
void
LAPACK_DECL(dorg2r)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dorgbr --------------------------------------------------------------------
void
LAPACK_DECL(dorgbr)(const char       *VECT,
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
LAPACK_DECL(dorghr)(const INTEGER    *N,
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
LAPACK_DECL(dorgl2)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dorglq --------------------------------------------------------------------
void
LAPACK_DECL(dorglq)(const INTEGER    *M,
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
LAPACK_DECL(dorgql)(const INTEGER    *M,
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
LAPACK_DECL(dorgqr)(const INTEGER    *M,
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
LAPACK_DECL(dorgr2)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dorgrq --------------------------------------------------------------------
void
LAPACK_DECL(dorgrq)(const INTEGER    *M,
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
LAPACK_DECL(dorgtr)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dorm2l --------------------------------------------------------------------
void
LAPACK_DECL(dorm2l)(const char       *SIDE,
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
LAPACK_DECL(dorm2r)(const char       *SIDE,
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
LAPACK_DECL(dormbr)(const char       *VECT,
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
LAPACK_DECL(dormhr)(const char       *SIDE,
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
LAPACK_DECL(dorml2)(const char       *SIDE,
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
LAPACK_DECL(dormlq)(const char       *SIDE,
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
LAPACK_DECL(dormql)(const char       *SIDE,
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
LAPACK_DECL(dormqr)(const char       *SIDE,
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
LAPACK_DECL(dormr2)(const char       *SIDE,
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
LAPACK_DECL(dormr3)(const char       *SIDE,
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
LAPACK_DECL(dormrq)(const char       *SIDE,
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
LAPACK_DECL(dormrz)(const char       *SIDE,
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
LAPACK_DECL(dormtr)(const char       *SIDE,
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
LAPACK_DECL(dpbcon)(const char       *UPLO,
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
LAPACK_DECL(dpbequ)(const char       *UPLO,
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
LAPACK_DECL(dpbrfs)(const char       *UPLO,
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
LAPACK_DECL(dpbstf)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO);

//-- dpbsv ---------------------------------------------------------------------
void
LAPACK_DECL(dpbsv)(const char           *UPLO,
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
LAPACK_DECL(dpbsvx)(const char       *FACT,
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
LAPACK_DECL(dpbtf2)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO);

//-- dpbtrf --------------------------------------------------------------------
void
LAPACK_DECL(dpbtrf)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO);

//-- dpbtrs --------------------------------------------------------------------
void
LAPACK_DECL(dpbtrs)(const char       *UPLO,
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
LAPACK_DECL(dpftrf)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    INTEGER          *INFO);

//-- dpftri --------------------------------------------------------------------
void
LAPACK_DECL(dpftri)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    INTEGER          *INFO);

//-- dpftrs --------------------------------------------------------------------
void
LAPACK_DECL(dpftrs)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dpocon --------------------------------------------------------------------
void
LAPACK_DECL(dpocon)(const char       *UPLO,
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
LAPACK_DECL(dpoequ)(const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *S,
                    DOUBLE           *SCOND,
                    DOUBLE           *AMAX,
                    INTEGER          *INFO);

//-- dpoequb -------------------------------------------------------------------
void
LAPACK_DECL(dpoequb)(const INTEGER    *N,
                     const DOUBLE     *A,
                     const INTEGER    *LDA,
                     DOUBLE           *S,
                     DOUBLE           *SCOND,
                     DOUBLE           *AMAX,
                     INTEGER          *INFO);

//-- dporfs --------------------------------------------------------------------
void
LAPACK_DECL(dporfs)(const char       *UPLO,
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
LAPACK_DECL(dposv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- dposvx --------------------------------------------------------------------
void
LAPACK_DECL(dposvx)(const char       *FACT,
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
LAPACK_DECL(dpotf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dpotrf --------------------------------------------------------------------
void
LAPACK_DECL(dpotrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dpotri --------------------------------------------------------------------
void
LAPACK_DECL(dpotri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dpotrs --------------------------------------------------------------------
void
LAPACK_DECL(dpotrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dppcon --------------------------------------------------------------------
void
LAPACK_DECL(dppcon)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO);

//-- dppequ --------------------------------------------------------------------
void
LAPACK_DECL(dppequ)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *S,
                    DOUBLE           *SCOND,
                    DOUBLE           *AMAX,
                    INTEGER          *INFO);

//-- dpprfs --------------------------------------------------------------------
void
LAPACK_DECL(dpprfs)(const char       *UPLO,
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
LAPACK_DECL(dppsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *AP,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- dppsvx --------------------------------------------------------------------
void
LAPACK_DECL(dppsvx)(const char       *FACT,
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
LAPACK_DECL(dpptrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    INTEGER          *INFO);

//-- dpptri --------------------------------------------------------------------
void
LAPACK_DECL(dpptri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    INTEGER          *INFO);

//-- dpptrs --------------------------------------------------------------------
void
LAPACK_DECL(dpptrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AP,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dpstf2 --------------------------------------------------------------------
void
LAPACK_DECL(dpstf2)(const char       *UPLO,
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
LAPACK_DECL(dpstrf)(const char       *UPLO,
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
LAPACK_DECL(dptcon)(const INTEGER    *N,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dpteqr --------------------------------------------------------------------
void
LAPACK_DECL(dpteqr)(const char       *COMPZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dptrfs --------------------------------------------------------------------
void
LAPACK_DECL(dptrfs)(const INTEGER    *N,
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
LAPACK_DECL(dptsv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *D,
                   DOUBLE               *E,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- dptsvx --------------------------------------------------------------------
void
LAPACK_DECL(dptsvx)(const char       *FACT,
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
LAPACK_DECL(dpttrf)(const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    INTEGER          *INFO);

//-- dpttrs --------------------------------------------------------------------
void
LAPACK_DECL(dpttrs)(const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dptts2 --------------------------------------------------------------------
void
LAPACK_DECL(dptts2)(const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    DOUBLE           *B,
                    const INTEGER    *LDB);

//-- drscl ---------------------------------------------------------------------
void
LAPACK_DECL(drscl)(const INTEGER        *N,
                   const DOUBLE         *SA,
                   DOUBLE               *SX,
                   const INTEGER        *INCX);

//-- dsbev ---------------------------------------------------------------------
void
LAPACK_DECL(dsbev)(const char           *JOBZ,
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
LAPACK_DECL(dsbevd)(const char       *JOBZ,
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
LAPACK_DECL(dsbevx)(const char       *JOBZ,
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
LAPACK_DECL(dsbgst)(const char       *VECT,
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
LAPACK_DECL(dsbgv)(const char           *JOBZ,
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
LAPACK_DECL(dsbgvd)(const char       *JOBZ,
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
LAPACK_DECL(dsbgvx)(const char       *JOBZ,
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
LAPACK_DECL(dsbtrd)(const char       *VECT,
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
LAPACK_DECL(dsecnd)();

//-- dsfrk ---------------------------------------------------------------------
void
LAPACK_DECL(dsfrk)(const char           *TRANSR,
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
LAPACK_DECL(dsgesv)(const INTEGER    *N,
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
LAPACK_DECL(dspcon)(const char       *UPLO,
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
LAPACK_DECL(dspev)(const char           *JOBZ,
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
LAPACK_DECL(dspevd)(const char       *JOBZ,
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
LAPACK_DECL(dspevx)(const char       *JOBZ,
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
LAPACK_DECL(dspgst)(const INTEGER    *ITYPE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    const DOUBLE     *BP,
                    INTEGER          *INFO);

//-- dspgv ---------------------------------------------------------------------
void
LAPACK_DECL(dspgv)(const INTEGER        *ITYPE,
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
LAPACK_DECL(dspgvd)(const INTEGER    *ITYPE,
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
LAPACK_DECL(dspgvx)(const INTEGER    *ITYPE,
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
LAPACK_DECL(dsposv)(const char       *UPLO,
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
LAPACK_DECL(dsprfs)(const char       *UPLO,
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
LAPACK_DECL(dspsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *AP,
                   INTEGER              *IPIV,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- dspsvx --------------------------------------------------------------------
void
LAPACK_DECL(dspsvx)(const char       *FACT,
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
LAPACK_DECL(dsptrd)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *TAU,
                    INTEGER          *INFO);

//-- dsptrf --------------------------------------------------------------------
void
LAPACK_DECL(dsptrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- dsptri --------------------------------------------------------------------
void
LAPACK_DECL(dsptri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    const INTEGER    *IPIV,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dsptrs --------------------------------------------------------------------
void
LAPACK_DECL(dsptrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AP,
                    const INTEGER    *IPIV,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dstebz --------------------------------------------------------------------
void
LAPACK_DECL(dstebz)(const char       *RANGE,
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
LAPACK_DECL(dstedc)(const char       *COMPZ,
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
LAPACK_DECL(dstegr)(const char       *JOBZ,
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
LAPACK_DECL(dstein)(const INTEGER    *N,
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
LAPACK_DECL(dstemr)(const char       *JOBZ,
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
LAPACK_DECL(dsteqr)(const char       *COMPZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dsterf --------------------------------------------------------------------
void
LAPACK_DECL(dsterf)(const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    INTEGER          *INFO);

//-- dstev ---------------------------------------------------------------------
void
LAPACK_DECL(dstev)(const char           *JOBZ,
                   const INTEGER        *N,
                   DOUBLE               *D,
                   DOUBLE               *E,
                   DOUBLE               *Z,
                   const INTEGER        *LDZ,
                   DOUBLE               *WORK,
                   INTEGER              *INFO);

//-- dstevd --------------------------------------------------------------------
void
LAPACK_DECL(dstevd)(const char       *JOBZ,
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
LAPACK_DECL(dstevr)(const char       *JOBZ,
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
LAPACK_DECL(dstevx)(const char       *JOBZ,
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
LAPACK_DECL(dsycon)(const char       *UPLO,
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
LAPACK_DECL(dsyconv)(const char       *UPLO,
                     const char       *WAY,
                     const INTEGER    *N,
                     const DOUBLE     *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE           *WORK,
                     INTEGER          *INFO);

//-- dsyequb -------------------------------------------------------------------
void
LAPACK_DECL(dsyequb)(const char       *UPLO,
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
LAPACK_DECL(dsyev)(const char           *JOBZ,
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
LAPACK_DECL(dsyevd)(const char       *JOBZ,
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
LAPACK_DECL(dsyevr)(const char       *JOBZ,
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
LAPACK_DECL(dsyevx)(const char       *JOBZ,
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
LAPACK_DECL(dsygs2)(const INTEGER    *ITYPE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dsygst --------------------------------------------------------------------
void
LAPACK_DECL(dsygst)(const INTEGER    *ITYPE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- dsygv ---------------------------------------------------------------------
void
LAPACK_DECL(dsygv)(const INTEGER        *ITYPE,
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
LAPACK_DECL(dsygvd)(const INTEGER    *ITYPE,
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
LAPACK_DECL(dsygvx)(const INTEGER    *ITYPE,
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
LAPACK_DECL(dsyrfs)(const char       *UPLO,
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
LAPACK_DECL(dsysv)(const char           *UPLO,
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
LAPACK_DECL(dsysvx)(const char       *FACT,
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
LAPACK_DECL(dsyswapr)(const char       *UPLO,
                      const INTEGER    *N,
                      DOUBLE           *A,
                      const INTEGER    *LDA,
                      const INTEGER    *I1,
                      const INTEGER    *I2);

//-- dsytd2 --------------------------------------------------------------------
void
LAPACK_DECL(dsytd2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *TAU,
                    INTEGER          *INFO);

//-- dsytf2 --------------------------------------------------------------------
void
LAPACK_DECL(dsytf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- dsytrd --------------------------------------------------------------------
void
LAPACK_DECL(dsytrd)(const char       *UPLO,
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
LAPACK_DECL(dsytrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dsytri --------------------------------------------------------------------
void
LAPACK_DECL(dsytri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- dsytri2 -------------------------------------------------------------------
void
LAPACK_DECL(dsytri2)(const char       *UPLO,
                     const INTEGER    *N,
                     DOUBLE           *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE           *WORK,
                     const INTEGER    *LWORK,
                     INTEGER          *INFO);

//-- dsytri2x ------------------------------------------------------------------
void
LAPACK_DECL(dsytri2x)(const char       *UPLO,
                      const INTEGER    *N,
                      DOUBLE           *A,
                      const INTEGER    *LDA,
                      const INTEGER    *IPIV,
                      DOUBLE           *WORK,
                      const INTEGER    *NB,
                      INTEGER          *INFO);

//-- dsytrs --------------------------------------------------------------------
void
LAPACK_DECL(dsytrs)(const char       *UPLO,
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
LAPACK_DECL(dsytrs2)(const char       *UPLO,
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
LAPACK_DECL(dtbcon)(const char       *NORM,
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
LAPACK_DECL(dtbrfs)(const char       *UPLO,
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
LAPACK_DECL(dtbtrs)(const char       *UPLO,
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
LAPACK_DECL(dtfsm)(const char           *TRANSR,
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
LAPACK_DECL(dtftri)(const char       *TRANSR,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    INTEGER          *INFO);

//-- dtfttp --------------------------------------------------------------------
void
LAPACK_DECL(dtfttp)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *ARF,
                    DOUBLE           *AP,
                    INTEGER          *INFO);

//-- dtfttr --------------------------------------------------------------------
void
LAPACK_DECL(dtfttr)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *ARF,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dtgevc --------------------------------------------------------------------
void
LAPACK_DECL(dtgevc)(const char       *SIDE,
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
LAPACK_DECL(dtgex2)(const LOGICAL    *WANTQ,
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
LAPACK_DECL(dtgexc)(const LOGICAL    *WANTQ,
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
LAPACK_DECL(dtgsen)(const INTEGER    *IJOB,
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
LAPACK_DECL(dtgsja)(const char       *JOBU,
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
LAPACK_DECL(dtgsna)(const char       *JOB,
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
LAPACK_DECL(dtgsy2)(const char       *TRANS,
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
LAPACK_DECL(dtgsyl)(const char       *TRANS,
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
LAPACK_DECL(dtpcon)(const char       *NORM,
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
LAPACK_DECL(dtprfs)(const char       *UPLO,
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
LAPACK_DECL(dtptri)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    INTEGER          *INFO);

//-- dtptrs --------------------------------------------------------------------
void
LAPACK_DECL(dtptrs)(const char       *UPLO,
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
LAPACK_DECL(dtpttf)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *ARF,
                    INTEGER          *INFO);

//-- dtpttr --------------------------------------------------------------------
void
LAPACK_DECL(dtpttr)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dtrcon --------------------------------------------------------------------
void
LAPACK_DECL(dtrcon)(const char       *NORM,
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
LAPACK_DECL(dtrevc)(const char       *SIDE,
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
LAPACK_DECL(dtrexc)(const char       *COMPQ,
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
LAPACK_DECL(dtrrfs)(const char       *UPLO,
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
LAPACK_DECL(dtrsen)(const char       *JOB,
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
LAPACK_DECL(dtrsna)(const char       *JOB,
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
LAPACK_DECL(dtrsyl)(const char       *TRANA,
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
LAPACK_DECL(dtrti2)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dtrtri --------------------------------------------------------------------
void
LAPACK_DECL(dtrtri)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- dtrtrs --------------------------------------------------------------------
void
LAPACK_DECL(dtrtrs)(const char       *UPLO,
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
LAPACK_DECL(dtrttf)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *ARF,
                    INTEGER          *INFO);

//-- dtrttp --------------------------------------------------------------------
void
LAPACK_DECL(dtrttp)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *AP,
                    INTEGER          *INFO);

//-- dtzrqf --------------------------------------------------------------------
void
LAPACK_DECL(dtzrqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    INTEGER          *INFO);

//-- dtzrzf --------------------------------------------------------------------
void
LAPACK_DECL(dtzrzf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- dzsum1 --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dzsum1)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *CX,
                    const INTEGER            *INCX);

//-- ieeeck --------------------------------------------------------------------
INTEGER
LAPACK_DECL(ieeeck)(const INTEGER    *ISPEC,
                    const FLOAT      *ZERO,
                    const FLOAT      *ONE);

//-- iladlc --------------------------------------------------------------------
INTEGER
LAPACK_DECL(iladlc)(const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA);

//-- iladlr --------------------------------------------------------------------
INTEGER
LAPACK_DECL(iladlr)(const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA);

//-- ilaslc --------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilaslc)(const INTEGER    *M,
                    const INTEGER    *N,
                    const FLOAT      *A,
                    const INTEGER    *LDA);

//-- ilaslr --------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilaslr)(const INTEGER    *M,
                    const INTEGER    *N,
                    const FLOAT      *A,
                    const INTEGER    *LDA);

//-- ilatrans ------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilatrans)(const char       *TRANS);

//-- ilauplo -------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilauplo)(const char   *UPLO);

//-- ilaver --------------------------------------------------------------------
void
LAPACK_DECL(ilaver)(INTEGER  *VERS_MAJOR,
                    INTEGER  *VERS_MINOR,
                    INTEGER  *VERS_PATCH);

//-- lsame ---------------------------------------------------------------------
LOGICAL
LAPACK_DECL(lsame)(const char       *CA,
                   const char       *CB);

//-- lsamen --------------------------------------------------------------------
LOGICAL
LAPACK_DECL(lsamen)(const INTEGER    *N,
                    const char       *CA,
                    const char       *CB);

//-- sgetf2 --------------------------------------------------------------------
void
LAPACK_DECL(sgetf2)(const INTEGER    *M,
                    const INTEGER    *N,
                    FLOAT            *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- sgetrf --------------------------------------------------------------------
void
LAPACK_DECL(sgetrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    FLOAT            *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- sgetrs --------------------------------------------------------------------
void
LAPACK_DECL(sgetrs)(const char       *TRANS,
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
LAPACK_DECL(sisnan)(const FLOAT  *SIN);

//-- slag2d --------------------------------------------------------------------
void
LAPACK_DECL(slag2d)(const INTEGER    *M,
                    const INTEGER    *N,
                    const FLOAT      *SA,
                    const INTEGER    *LDSA,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- slaisnan ------------------------------------------------------------------
LOGICAL
LAPACK_DECL(slaisnan)(const FLOAT      *SIN1,
                      const FLOAT      *SIN2);

//-- slamch --------------------------------------------------------------------
FLOAT
LAPACK_DECL(slamch)(const char   *CMACH);

//-- slaswp --------------------------------------------------------------------
void
LAPACK_DECL(slaswp)(const INTEGER    *N,
                    FLOAT            *A,
                    const INTEGER    *LDA,
                    const INTEGER    *K1,
                    const INTEGER    *K2,
                    const INTEGER    *IPIV,
                    const INTEGER    *INCX);

//-- spotf2 --------------------------------------------------------------------
void
LAPACK_DECL(spotf2)(const char       *UPLO,
                    const INTEGER    *N,
                    FLOAT            *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- spotrf --------------------------------------------------------------------
void
LAPACK_DECL(spotrf)(const char       *UPLO,
                    const INTEGER    *N,
                    FLOAT            *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- spotrs --------------------------------------------------------------------
void
LAPACK_DECL(spotrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const FLOAT      *A,
                    const INTEGER    *LDA,
                    FLOAT            *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- zbbcsd --------------------------------------------------------------------
void
LAPACK_DECL(zbbcsd)(const char       *JOBU1,
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
LAPACK_DECL(zbdsqr)(const char       *UPLO,
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
LAPACK_DECL(zcgesv)(const INTEGER            *N,
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
LAPACK_DECL(zcposv)(const char               *UPLO,
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
LAPACK_DECL(zdrscl)(const INTEGER    *N,
                    const DOUBLE     *SA,
                    DOUBLE_COMPLEX   *SX,
                    const INTEGER    *INCX);

//-- zgbbrd --------------------------------------------------------------------
void
LAPACK_DECL(zgbbrd)(const char       *VECT,
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
LAPACK_DECL(zgbcon)(const char               *NORM,
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
LAPACK_DECL(zgbequ)(const INTEGER            *M,
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
LAPACK_DECL(zgbequb)(const INTEGER            *M,
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
LAPACK_DECL(zgbrfs)(const char               *TRANS,
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
LAPACK_DECL(zgbrfsx)(const char               *TRANS,
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
LAPACK_DECL(zgbsv)(const INTEGER        *N,
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
LAPACK_DECL(zgbsvx)(const char       *FACT,
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
LAPACK_DECL(zgbsvxx)(const char       *FACT,
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
LAPACK_DECL(zgbtf2)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zgbtrf --------------------------------------------------------------------
void
LAPACK_DECL(zgbtrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zgbtrs --------------------------------------------------------------------
void
LAPACK_DECL(zgbtrs)(const char               *TRANS,
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
LAPACK_DECL(zgebak)(const char       *JOB,
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
LAPACK_DECL(zgebal)(const char       *JOB,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *ILO,
                    INTEGER          *IHI,
                    DOUBLE           *SCALE,
                    INTEGER          *INFO);

//-- zgebd2 --------------------------------------------------------------------
void
LAPACK_DECL(zgebd2)(const INTEGER    *M,
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
LAPACK_DECL(zgebrd)(const INTEGER    *M,
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
LAPACK_DECL(zgecon)(const char               *NORM,
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
LAPACK_DECL(zgeequ)(const INTEGER            *M,
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
LAPACK_DECL(zgeequb)(const INTEGER            *M,
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
LAPACK_DECL(zgees)(const char           *JOBVS,
                   const char           *SORT,
                   const LOGICAL        *SELECT,
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
LAPACK_DECL(zgeesx)(const char       *JOBVS,
                    const char       *SORT,
                    const LOGICAL    *SELECT,
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
LAPACK_DECL(zgeev)(const char           *JOBVL,
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
LAPACK_DECL(zgeevx)(const char       *BALANC,
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
LAPACK_DECL(zgegs)(const char           *JOBVSL,
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
LAPACK_DECL(zgegv)(const char           *JOBVL,
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
LAPACK_DECL(zgehd2)(const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zgehrd --------------------------------------------------------------------
void
LAPACK_DECL(zgehrd)(const INTEGER    *N,
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
LAPACK_DECL(zgelq2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zgelqf --------------------------------------------------------------------
void
LAPACK_DECL(zgelqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zgels ---------------------------------------------------------------------
void
LAPACK_DECL(zgels)(const char           *TRANS,
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
LAPACK_DECL(zgelsd)(const INTEGER            *M,
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
LAPACK_DECL(zgelss)(const INTEGER    *M,
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
LAPACK_DECL(zgelsx)(const INTEGER    *M,
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
LAPACK_DECL(zgelsy)(const INTEGER    *M,
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
LAPACK_DECL(zgeql2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zgeqlf --------------------------------------------------------------------
void
LAPACK_DECL(zgeqlf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zgeqp3 --------------------------------------------------------------------
void
LAPACK_DECL(zgeqp3)(const INTEGER    *M,
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
LAPACK_DECL(zgeqpf)(const INTEGER    *M,
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
LAPACK_DECL(zgeqr2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zgeqr2p -------------------------------------------------------------------
void
LAPACK_DECL(zgeqr2p)(const INTEGER    *M,
                     const INTEGER    *N,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     DOUBLE_COMPLEX   *TAU,
                     DOUBLE_COMPLEX   *WORK,
                     INTEGER          *INFO);

//-- zgeqrf --------------------------------------------------------------------
void
LAPACK_DECL(zgeqrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zgeqrfp -------------------------------------------------------------------
void
LAPACK_DECL(zgeqrfp)(const INTEGER    *M,
                     const INTEGER    *N,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     DOUBLE_COMPLEX   *TAU,
                     DOUBLE_COMPLEX   *WORK,
                     const INTEGER    *LWORK,
                     INTEGER          *INFO);

//-- zgerfs --------------------------------------------------------------------
void
LAPACK_DECL(zgerfs)(const char               *TRANS,
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
LAPACK_DECL(zgerfsx)(const char               *TRANS,
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
LAPACK_DECL(zgerq2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zgerqf --------------------------------------------------------------------
void
LAPACK_DECL(zgerqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zgesc2 --------------------------------------------------------------------
void
LAPACK_DECL(zgesc2)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *RHS,
                    const INTEGER            *IPIV,
                    const INTEGER            *JPIV,
                    DOUBLE                   *SCALE);

//-- zgesdd --------------------------------------------------------------------
void
LAPACK_DECL(zgesdd)(const char       *JOBZ,
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
LAPACK_DECL(zgesv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   INTEGER              *IPIV,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zgesvd --------------------------------------------------------------------
void
LAPACK_DECL(zgesvd)(const char       *JOBU,
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
LAPACK_DECL(zgesvx)(const char       *FACT,
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
LAPACK_DECL(zgesvxx)(const char       *FACT,
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
LAPACK_DECL(zgetc2)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *JPIV,
                    INTEGER          *INFO);

//-- zgetf2 --------------------------------------------------------------------
void
LAPACK_DECL(zgetf2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zgetrf --------------------------------------------------------------------
void
LAPACK_DECL(zgetrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zgetri --------------------------------------------------------------------
void
LAPACK_DECL(zgetri)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zgetrs --------------------------------------------------------------------
void
LAPACK_DECL(zgetrs)(const char               *TRANS,
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
LAPACK_DECL(zggbak)(const char       *JOB,
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
LAPACK_DECL(zggbal)(const char       *JOB,
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
LAPACK_DECL(zgges)(const char           *JOBVSL,
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
LAPACK_DECL(zggesx)(const char       *JOBVSL,
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
LAPACK_DECL(zggev)(const char           *JOBVL,
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
LAPACK_DECL(zggevx)(const char       *BALANC,
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
LAPACK_DECL(zggglm)(const INTEGER    *N,
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
LAPACK_DECL(zgghrd)(const char       *COMPQ,
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
LAPACK_DECL(zgglse)(const INTEGER    *M,
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
LAPACK_DECL(zggqrf)(const INTEGER    *N,
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
LAPACK_DECL(zggrqf)(const INTEGER    *M,
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
LAPACK_DECL(zggsvd)(const char       *JOBU,
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
LAPACK_DECL(zggsvp)(const char       *JOBU,
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
LAPACK_DECL(zgtcon)(const char               *NORM,
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
LAPACK_DECL(zgtrfs)(const char               *TRANS,
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
LAPACK_DECL(zgtsv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *DL,
                   DOUBLE_COMPLEX       *D,
                   DOUBLE_COMPLEX       *DU,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zgtsvx --------------------------------------------------------------------
void
LAPACK_DECL(zgtsvx)(const char               *FACT,
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
LAPACK_DECL(zgttrf)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *DL,
                    DOUBLE_COMPLEX   *D,
                    DOUBLE_COMPLEX   *DU,
                    DOUBLE_COMPLEX   *DU2,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zgttrs --------------------------------------------------------------------
void
LAPACK_DECL(zgttrs)(const char               *TRANS,
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
LAPACK_DECL(zgtts2)(const INTEGER            *ITRANS,
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
LAPACK_DECL(zhbev)(const char           *JOBZ,
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
LAPACK_DECL(zhbevd)(const char       *JOBZ,
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
LAPACK_DECL(zhbevx)(const char       *JOBZ,
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
LAPACK_DECL(zhbgst)(const char               *VECT,
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
LAPACK_DECL(zhbgv)(const char           *JOBZ,
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
LAPACK_DECL(zhbgvd)(const char       *JOBZ,
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
LAPACK_DECL(zhbgvx)(const char       *JOBZ,
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
LAPACK_DECL(zhbtrd)(const char       *VECT,
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
LAPACK_DECL(zhecon)(const char               *UPLO,
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
LAPACK_DECL(zheequb)(const char               *UPLO,
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
LAPACK_DECL(zheev)(const char           *JOBZ,
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
LAPACK_DECL(zheevd)(const char       *JOBZ,
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
LAPACK_DECL(zheevr)(const char       *JOBZ,
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
LAPACK_DECL(zheevx)(const char       *JOBZ,
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
LAPACK_DECL(zhegs2)(const INTEGER            *ITYPE,
                    const char               *UPLO,
                    const INTEGER            *N,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zhegst --------------------------------------------------------------------
void
LAPACK_DECL(zhegst)(const INTEGER            *ITYPE,
                    const char               *UPLO,
                    const INTEGER            *N,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zhegv ---------------------------------------------------------------------
void
LAPACK_DECL(zhegv)(const INTEGER        *ITYPE,
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
LAPACK_DECL(zhegvd)(const INTEGER    *ITYPE,
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
LAPACK_DECL(zhegvx)(const INTEGER    *ITYPE,
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
LAPACK_DECL(zherfs)(const char               *UPLO,
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
LAPACK_DECL(zherfsx)(const char               *UPLO,
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
LAPACK_DECL(zhesv)(const char           *UPLO,
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
LAPACK_DECL(zhesvx)(const char               *FACT,
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
LAPACK_DECL(zhesvxx)(const char       *FACT,
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
LAPACK_DECL(zheswapr)(const char           *UPLO,
                      const INTEGER        *N,
                      DOUBLE_COMPLEX       *A,
                      const INTEGER        *LDA,
                      const INTEGER        *I1,
                      const INTEGER        *I2);

//-- zhetd2 --------------------------------------------------------------------
void
LAPACK_DECL(zhetd2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *TAU,
                    INTEGER          *INFO);

//-- zhetf2 --------------------------------------------------------------------
void
LAPACK_DECL(zhetf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zhetrd --------------------------------------------------------------------
void
LAPACK_DECL(zhetrd)(const char       *UPLO,
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
LAPACK_DECL(zhetrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zhetri --------------------------------------------------------------------
void
LAPACK_DECL(zhetri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zhetri2 -------------------------------------------------------------------
void
LAPACK_DECL(zhetri2)(const char       *UPLO,
                     const INTEGER    *N,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE_COMPLEX   *WORK,
                     const INTEGER    *LWORK,
                     INTEGER          *INFO);

//-- zhetri2x ------------------------------------------------------------------
void
LAPACK_DECL(zhetri2x)(const char           *UPLO,
                      const INTEGER        *N,
                      DOUBLE_COMPLEX       *A,
                      const INTEGER        *LDA,
                      const INTEGER        *IPIV,
                      DOUBLE_COMPLEX       *WORK,
                      const INTEGER        *NB,
                      INTEGER              *INFO);

//-- zhetrs --------------------------------------------------------------------
void
LAPACK_DECL(zhetrs)(const char               *UPLO,
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
LAPACK_DECL(zhetrs2)(const char       *UPLO,
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
LAPACK_DECL(zhfrk)(const char               *TRANSR,
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
LAPACK_DECL(zhgeqz)(const char       *JOB,
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
LAPACK_DECL(zhpcon)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    const INTEGER            *IPIV,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zhpev ---------------------------------------------------------------------
void
LAPACK_DECL(zhpev)(const char           *JOBZ,
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
LAPACK_DECL(zhpevd)(const char       *JOBZ,
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
LAPACK_DECL(zhpevx)(const char       *JOBZ,
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
LAPACK_DECL(zhpgst)(const INTEGER            *ITYPE,
                    const char               *UPLO,
                    const INTEGER            *N,
                    DOUBLE_COMPLEX           *AP,
                    const DOUBLE_COMPLEX     *BP,
                    INTEGER                  *INFO);

//-- zhpgv ---------------------------------------------------------------------
void
LAPACK_DECL(zhpgv)(const INTEGER        *ITYPE,
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
LAPACK_DECL(zhpgvd)(const INTEGER    *ITYPE,
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
LAPACK_DECL(zhpgvx)(const INTEGER    *ITYPE,
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
LAPACK_DECL(zhprfs)(const char               *UPLO,
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
LAPACK_DECL(zhpsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *AP,
                   INTEGER              *IPIV,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zhpsvx --------------------------------------------------------------------
void
LAPACK_DECL(zhpsvx)(const char               *FACT,
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
LAPACK_DECL(zhptrd)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *TAU,
                    INTEGER          *INFO);

//-- zhptrf --------------------------------------------------------------------
void
LAPACK_DECL(zhptrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zhptri --------------------------------------------------------------------
void
LAPACK_DECL(zhptri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    const INTEGER    *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zhptrs --------------------------------------------------------------------
void
LAPACK_DECL(zhptrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    const INTEGER            *IPIV,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zhsein --------------------------------------------------------------------
void
LAPACK_DECL(zhsein)(const char               *SIDE,
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
LAPACK_DECL(zhseqr)(const char       *JOB,
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
LAPACK_DECL(zla_gbamv)(const INTEGER            *TRANS,
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
LAPACK_DECL(zla_gbrcond_c)(const char               *TRANS,
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
LAPACK_DECL(zla_gbrcond_x)(const char               *TRANS,
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
LAPACK_DECL(zla_gbrfsx_extended)(const INTEGER            *PREC_TYPE,
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
LAPACK_DECL(zla_gbrpvgrw)(const INTEGER            *N,
                          const INTEGER            *KL,
                          const INTEGER            *KU,
                          const INTEGER            *NCOLS,
                          const DOUBLE_COMPLEX     *AB,
                          const INTEGER            *LDAB,
                          const DOUBLE_COMPLEX     *AFB,
                          const INTEGER            *LDAFB);

//-- zla_geamv -----------------------------------------------------------------
void
LAPACK_DECL(zla_geamv)(const INTEGER            *TRANS,
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
LAPACK_DECL(zla_gercond_c)(const char               *TRANS,
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
LAPACK_DECL(zla_gercond_x)(const char               *TRANS,
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
LAPACK_DECL(zla_gerfsx_extended)(const INTEGER            *PREC_TYPE,
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
LAPACK_DECL(zla_heamv)(const INTEGER            *UPLO,
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
LAPACK_DECL(zla_hercond_c)(const char               *UPLO,
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
LAPACK_DECL(zla_hercond_x)(const char               *UPLO,
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
LAPACK_DECL(zla_herfsx_extended)(const INTEGER            *PREC_TYPE,
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
LAPACK_DECL(zla_herpvgrw)(const char               *UPLO,
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
LAPACK_DECL(zla_lin_berr)(const INTEGER            *N,
                          const INTEGER            *NZ,
                          const INTEGER            *NRHS,
                          const DOUBLE_COMPLEX     *RES,
                          const DOUBLE             *AYB,
                          DOUBLE                   *BERR);

//-- zla_porcond_c -------------------------------------------------------------
DOUBLE
LAPACK_DECL(zla_porcond_c)(const char               *UPLO,
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
LAPACK_DECL(zla_porcond_x)(const char               *UPLO,
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
LAPACK_DECL(zla_porfsx_extended)(const INTEGER            *PREC_TYPE,
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
LAPACK_DECL(zla_porpvgrw)(const char               *UPLO,
                          const INTEGER            *NCOLS,
                          const DOUBLE_COMPLEX     *A,
                          const INTEGER            *LDA,
                          const DOUBLE_COMPLEX     *AF,
                          const INTEGER            *LDAF,
                          const DOUBLE             *WORK);

//-- zla_rpvgrw ----------------------------------------------------------------
DOUBLE
LAPACK_DECL(zla_rpvgrw)(const INTEGER            *N,
                        const INTEGER            *NCOLS,
                        const DOUBLE_COMPLEX     *A,
                        const INTEGER            *LDA,
                        const DOUBLE_COMPLEX     *AF,
                        const INTEGER            *LDAF);

//-- zla_syamv -----------------------------------------------------------------
void
LAPACK_DECL(zla_syamv)(const INTEGER            *UPLO,
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
LAPACK_DECL(zla_syrcond_c)(const char               *UPLO,
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
LAPACK_DECL(zla_syrcond_x)(const char               *UPLO,
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
LAPACK_DECL(zla_syrfsx_extended)(const INTEGER            *PREC_TYPE,
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
LAPACK_DECL(zla_syrpvgrw)(const char               *UPLO,
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
LAPACK_DECL(zla_wwaddw)(const INTEGER            *N,
                        DOUBLE_COMPLEX           *X,
                        DOUBLE_COMPLEX           *Y,
                        const DOUBLE_COMPLEX     *W);

//-- zlabrd --------------------------------------------------------------------
void
LAPACK_DECL(zlabrd)(const INTEGER    *M,
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
LAPACK_DECL(zlacgv)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *INCX);

//-- zlacn2 --------------------------------------------------------------------
void
LAPACK_DECL(zlacn2)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *V,
                    DOUBLE_COMPLEX   *X,
                    DOUBLE           *EST,
                    INTEGER          *KASE,
                    INTEGER          *ISAVE);

//-- zlacon --------------------------------------------------------------------
void
LAPACK_DECL(zlacon)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *V,
                    DOUBLE_COMPLEX   *X,
                    DOUBLE           *EST,
                    INTEGER          *KASE);

//-- zlacp2 --------------------------------------------------------------------
void
LAPACK_DECL(zlacp2)(const char       *UPLO,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB);

//-- zlacpy --------------------------------------------------------------------
void
LAPACK_DECL(zlacpy)(const char               *UPLO,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB);

//-- zlacrm --------------------------------------------------------------------
void
LAPACK_DECL(zlacrm)(const INTEGER            *M,
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
LAPACK_DECL(zlacrt)(const INTEGER            *N,
                    DOUBLE_COMPLEX           *CX,
                    const INTEGER            *INCX,
                    DOUBLE_COMPLEX           *CY,
                    const INTEGER            *INCY,
                    const DOUBLE_COMPLEX     *C,
                    const DOUBLE_COMPLEX     *S);

//-- zladiv --------------------------------------------------------------------
/*
void
LAPACK_DECL(zladiv)(DOUBLE_COMPLEX           *ret,
                    const DOUBLE_COMPLEX     *X,
                    const DOUBLE_COMPLEX     *Y);
*/

//-- zlaed0 --------------------------------------------------------------------
void
LAPACK_DECL(zlaed0)(const INTEGER    *QSIZ,
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
LAPACK_DECL(zlaed7)(const INTEGER    *N,
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
LAPACK_DECL(zlaed8)(INTEGER          *K,
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
LAPACK_DECL(zlaein)(const LOGICAL            *RIGHTV,
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
LAPACK_DECL(zlaesy)(const DOUBLE_COMPLEX     *A,
                    const DOUBLE_COMPLEX     *B,
                    const DOUBLE_COMPLEX     *C,
                    DOUBLE_COMPLEX           *RT1,
                    DOUBLE_COMPLEX           *RT2,
                    DOUBLE_COMPLEX           *EVSCAL,
                    DOUBLE_COMPLEX           *CS1,
                    DOUBLE_COMPLEX           *SN1);

//-- zlaev2 --------------------------------------------------------------------
void
LAPACK_DECL(zlaev2)(const DOUBLE_COMPLEX     *A,
                    const DOUBLE_COMPLEX     *B,
                    const DOUBLE_COMPLEX     *C,
                    DOUBLE                   *RT1,
                    DOUBLE                   *RT2,
                    DOUBLE                   *CS1,
                    DOUBLE_COMPLEX           *SN1);

//-- zlag2c --------------------------------------------------------------------
void
LAPACK_DECL(zlag2c)(const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    FLOAT_COMPLEX            *SA,
                    const INTEGER            *LDSA,
                    INTEGER                  *INFO);

//-- zlags2 --------------------------------------------------------------------
void
LAPACK_DECL(zlags2)(const LOGICAL            *UPPER,
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
LAPACK_DECL(zlagtm)(const char               *TRANS,
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
LAPACK_DECL(zlahef)(const char       *UPLO,
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
LAPACK_DECL(zlahqr)(const LOGICAL    *WANTT,
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
LAPACK_DECL(zlahr2)(const INTEGER    *N,
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
LAPACK_DECL(zlahrd)(const INTEGER    *N,
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
LAPACK_DECL(zlaic1)(const INTEGER            *JOB,
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
LAPACK_DECL(zlals0)(const INTEGER    *ICOMPQ,
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
LAPACK_DECL(zlalsa)(const INTEGER    *ICOMPQ,
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
LAPACK_DECL(zlalsd)(const char       *UPLO,
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
LAPACK_DECL(zlangb)(const char               *NORM,
                    const INTEGER            *N,
                    const INTEGER            *KL,
                    const INTEGER            *KU,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE                   *WORK);

//-- zlange --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlange)(const char               *NORM,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *WORK);

//-- zlangt --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlangt)(const char               *NORM,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *DL,
                    const DOUBLE_COMPLEX     *D,
                    const DOUBLE_COMPLEX     *DU);

//-- zlanhb --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlanhb)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE                   *WORK);

//-- zlanhe --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlanhe)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *WORK);

//-- zlanhf --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlanhf)(const char               *NORM,
                    const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    DOUBLE                   *WORK);

//-- zlanhp --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlanhp)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE                   *WORK);

//-- zlanhs --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlanhs)(const char               *NORM,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *WORK);

//-- zlanht --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlanht)(const char               *NORM,
                    const INTEGER            *N,
                    const DOUBLE             *D,
                    const DOUBLE_COMPLEX     *E);

//-- zlansb --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlansb)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE                   *WORK);

//-- zlansp --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlansp)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE                   *WORK);

//-- zlansy --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlansy)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *WORK);

//-- zlantb --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlantb)(const char               *NORM,
                    const char               *UPLO,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE                   *WORK);

//-- zlantp --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlantp)(const char               *NORM,
                    const char               *UPLO,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE                   *WORK);

//-- zlantr --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlantr)(const char               *NORM,
                    const char               *UPLO,
                    const char               *DIAG,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *WORK);

//-- zlapll --------------------------------------------------------------------
void
LAPACK_DECL(zlapll)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *INCX,
                    DOUBLE_COMPLEX   *Y,
                    const INTEGER    *INCY,
                    DOUBLE           *SSMIN);

//-- zlapmr --------------------------------------------------------------------
void
LAPACK_DECL(zlapmr)(const LOGICAL    *FORWRD,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *LDX,
                    INTEGER          *K);

//-- zlapmt --------------------------------------------------------------------
void
LAPACK_DECL(zlapmt)(const LOGICAL    *FORWRD,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *LDX,
                    INTEGER          *K);

//-- zlaqgb --------------------------------------------------------------------
void
LAPACK_DECL(zlaqgb)(const INTEGER    *M,
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
LAPACK_DECL(zlaqge)(const INTEGER    *M,
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
LAPACK_DECL(zlaqhb)(const char       *UPLO,
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
LAPACK_DECL(zlaqhe)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- zlaqhp --------------------------------------------------------------------
void
LAPACK_DECL(zlaqhp)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- zlaqp2 --------------------------------------------------------------------
void
LAPACK_DECL(zlaqp2)(const INTEGER    *M,
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
LAPACK_DECL(zlaqps)(const INTEGER    *M,
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
LAPACK_DECL(zlaqr0)(const LOGICAL    *WANTT,
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
LAPACK_DECL(zlaqr1)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *H,
                    const INTEGER            *LDH,
                    const DOUBLE_COMPLEX     *S1,
                    const DOUBLE_COMPLEX     *S2,
                    DOUBLE_COMPLEX           *V);

//-- zlaqr2 --------------------------------------------------------------------
void
LAPACK_DECL(zlaqr2)(const LOGICAL    *WANTT,
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
LAPACK_DECL(zlaqr3)(const LOGICAL    *WANTT,
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
LAPACK_DECL(zlaqr4)(const LOGICAL    *WANTT,
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
LAPACK_DECL(zlaqr5)(const LOGICAL    *WANTT,
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
LAPACK_DECL(zlaqsb)(const char       *UPLO,
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
LAPACK_DECL(zlaqsp)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- zlaqsy --------------------------------------------------------------------
void
LAPACK_DECL(zlaqsy)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED);

//-- zlar1v --------------------------------------------------------------------
void
LAPACK_DECL(zlar1v)(const INTEGER    *N,
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
LAPACK_DECL(zlar2v)(const INTEGER            *N,
                    DOUBLE_COMPLEX           *X,
                    DOUBLE_COMPLEX           *Y,
                    DOUBLE_COMPLEX           *Z,
                    const INTEGER            *INCX,
                    const DOUBLE             *C,
                    const DOUBLE_COMPLEX     *S,
                    const INTEGER            *INCC);

//-- zlarcm --------------------------------------------------------------------
void
LAPACK_DECL(zlarcm)(const INTEGER            *M,
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
LAPACK_DECL(zlarf)(const char               *SIDE,
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
LAPACK_DECL(zlarfb)(const char               *SIDE,
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
LAPACK_DECL(zlarfg)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *ALPHA,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *INCX,
                    DOUBLE_COMPLEX   *TAU);

//-- zlarfgp -------------------------------------------------------------------
void
LAPACK_DECL(zlarfgp)(const INTEGER    *N,
                     DOUBLE_COMPLEX   *ALPHA,
                     DOUBLE_COMPLEX   *X,
                     const INTEGER    *INCX,
                     DOUBLE_COMPLEX   *TAU);

//-- zlarft --------------------------------------------------------------------
void
LAPACK_DECL(zlarft)(const char               *DIRECT,
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
LAPACK_DECL(zlarfx)(const char               *SIDE,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *V,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK);

//-- zlargv --------------------------------------------------------------------
void
LAPACK_DECL(zlargv)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *INCX,
                    DOUBLE_COMPLEX   *Y,
                    const INTEGER    *INCY,
                    DOUBLE           *C,
                    const INTEGER    *INCC);

//-- zlarnv --------------------------------------------------------------------
void
LAPACK_DECL(zlarnv)(const INTEGER    *IDIST,
                    INTEGER          *ISEED,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *X);

//-- zlarrv --------------------------------------------------------------------
void
LAPACK_DECL(zlarrv)(const INTEGER    *N,
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
LAPACK_DECL(zlarscl2)(const INTEGER        *M,
                      const INTEGER        *N,
                      const DOUBLE         *D,
                      DOUBLE_COMPLEX       *X,
                      const INTEGER        *LDX);

//-- zlartg --------------------------------------------------------------------
void
LAPACK_DECL(zlartg)(const DOUBLE_COMPLEX     *F,
                    const DOUBLE_COMPLEX     *G,
                    DOUBLE                   *CS,
                    DOUBLE_COMPLEX           *SN,
                    DOUBLE_COMPLEX           *R);

//-- zlartv --------------------------------------------------------------------
void
LAPACK_DECL(zlartv)(const INTEGER            *N,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *INCX,
                    DOUBLE_COMPLEX           *Y,
                    const INTEGER            *INCY,
                    const DOUBLE             *C,
                    const DOUBLE_COMPLEX     *S,
                    const INTEGER            *INCC);

//-- zlarz ---------------------------------------------------------------------
void
LAPACK_DECL(zlarz)(const char               *SIDE,
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
LAPACK_DECL(zlarzb)(const char               *SIDE,
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
LAPACK_DECL(zlarzt)(const char               *DIRECT,
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
LAPACK_DECL(zlascl)(const char       *TYPE,
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
LAPACK_DECL(zlascl2)(const INTEGER    *M,
                     const INTEGER    *N,
                     const DOUBLE     *D,
                     DOUBLE_COMPLEX   *X,
                     const INTEGER    *LDX);

//-- zlaset --------------------------------------------------------------------
void
LAPACK_DECL(zlaset)(const char               *UPLO,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *ALPHA,
                    const DOUBLE_COMPLEX     *BETA,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA);

//-- zlasr ---------------------------------------------------------------------
void
LAPACK_DECL(zlasr)(const char           *SIDE,
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
LAPACK_DECL(zlassq)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *X,
                    const INTEGER            *INCX,
                    DOUBLE                   *SCALE,
                    DOUBLE                   *SUMSQ);

//-- zlaswp --------------------------------------------------------------------
void
LAPACK_DECL(zlaswp)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const INTEGER    *K1,
                    const INTEGER    *K2,
                    const INTEGER    *IPIV,
                    const INTEGER    *INCX);

//-- zlasyf --------------------------------------------------------------------
void
LAPACK_DECL(zlasyf)(const char       *UPLO,
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
LAPACK_DECL(zlat2c)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    FLOAT_COMPLEX            *SA,
                    const INTEGER            *LDSA,
                    INTEGER                  *INFO);

//-- zlatbs --------------------------------------------------------------------
void
LAPACK_DECL(zlatbs)(const char               *UPLO,
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
LAPACK_DECL(zlatdf)(const INTEGER            *IJOB,
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
LAPACK_DECL(zlatps)(const char               *UPLO,
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
LAPACK_DECL(zlatrd)(const char       *UPLO,
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
LAPACK_DECL(zlatrs)(const char               *UPLO,
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
LAPACK_DECL(zlatrz)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *L,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK);

//-- zlatzm --------------------------------------------------------------------
void
LAPACK_DECL(zlatzm)(const char               *SIDE,
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
LAPACK_DECL(zlauu2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- zlauum --------------------------------------------------------------------
void
LAPACK_DECL(zlauum)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- zpbcon --------------------------------------------------------------------
void
LAPACK_DECL(zpbcon)(const char               *UPLO,
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
LAPACK_DECL(zpbequ)(const char               *UPLO,
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
LAPACK_DECL(zpbrfs)(const char               *UPLO,
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
LAPACK_DECL(zpbstf)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO);

//-- zpbsv ---------------------------------------------------------------------
void
LAPACK_DECL(zpbsv)(const char           *UPLO,
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
LAPACK_DECL(zpbsvx)(const char       *FACT,
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
LAPACK_DECL(zpbtf2)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO);

//-- zpbtrf --------------------------------------------------------------------
void
LAPACK_DECL(zpbtrf)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO);

//-- zpbtrs --------------------------------------------------------------------
void
LAPACK_DECL(zpbtrs)(const char               *UPLO,
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
LAPACK_DECL(zpftrf)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    INTEGER          *INFO);

//-- zpftri --------------------------------------------------------------------
void
LAPACK_DECL(zpftri)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    INTEGER          *INFO);

//-- zpftrs --------------------------------------------------------------------
void
LAPACK_DECL(zpftrs)(const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zpocon --------------------------------------------------------------------
void
LAPACK_DECL(zpocon)(const char               *UPLO,
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
LAPACK_DECL(zpoequ)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *S,
                    DOUBLE                   *SCOND,
                    DOUBLE                   *AMAX,
                    INTEGER                  *INFO);

//-- zpoequb -------------------------------------------------------------------
void
LAPACK_DECL(zpoequb)(const INTEGER            *N,
                     const DOUBLE_COMPLEX     *A,
                     const INTEGER            *LDA,
                     DOUBLE                   *S,
                     DOUBLE                   *SCOND,
                     DOUBLE                   *AMAX,
                     INTEGER                  *INFO);

//-- zporfs --------------------------------------------------------------------
void
LAPACK_DECL(zporfs)(const char               *UPLO,
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
LAPACK_DECL(zporfsx)(const char               *UPLO,
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
LAPACK_DECL(zposv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zposvx --------------------------------------------------------------------
void
LAPACK_DECL(zposvx)(const char       *FACT,
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
LAPACK_DECL(zposvxx)(const char       *FACT,
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
LAPACK_DECL(zpotf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- zpotrf --------------------------------------------------------------------
void
LAPACK_DECL(zpotrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- zpotri --------------------------------------------------------------------
void
LAPACK_DECL(zpotri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- zpotrs --------------------------------------------------------------------
void
LAPACK_DECL(zpotrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zppcon --------------------------------------------------------------------
void
LAPACK_DECL(zppcon)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zppequ --------------------------------------------------------------------
void
LAPACK_DECL(zppequ)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE                   *S,
                    DOUBLE                   *SCOND,
                    DOUBLE                   *AMAX,
                    INTEGER                  *INFO);

//-- zpprfs --------------------------------------------------------------------
void
LAPACK_DECL(zpprfs)(const char               *UPLO,
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
LAPACK_DECL(zppsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *AP,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zppsvx --------------------------------------------------------------------
void
LAPACK_DECL(zppsvx)(const char       *FACT,
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
LAPACK_DECL(zpptrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    INTEGER          *INFO);

//-- zpptri --------------------------------------------------------------------
void
LAPACK_DECL(zpptri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    INTEGER          *INFO);

//-- zpptrs --------------------------------------------------------------------
void
LAPACK_DECL(zpptrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zpstf2 --------------------------------------------------------------------
void
LAPACK_DECL(zpstf2)(const char       *UPLO,
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
LAPACK_DECL(zpstrf)(const char       *UPLO,
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
LAPACK_DECL(zptcon)(const INTEGER            *N,
                    const DOUBLE             *D,
                    const DOUBLE_COMPLEX     *E,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO);

//-- zpteqr --------------------------------------------------------------------
void
LAPACK_DECL(zpteqr)(const char       *COMPZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- zptrfs --------------------------------------------------------------------
void
LAPACK_DECL(zptrfs)(const char               *UPLO,
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
LAPACK_DECL(zptsv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *D,
                   DOUBLE_COMPLEX       *E,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zptsvx --------------------------------------------------------------------
void
LAPACK_DECL(zptsvx)(const char               *FACT,
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
LAPACK_DECL(zpttrf)(const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE_COMPLEX   *E,
                    INTEGER          *INFO);

//-- zpttrs --------------------------------------------------------------------
void
LAPACK_DECL(zpttrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE             *D,
                    const DOUBLE_COMPLEX     *E,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zptts2 --------------------------------------------------------------------
void
LAPACK_DECL(zptts2)(const INTEGER            *IUPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE             *D,
                    const DOUBLE_COMPLEX     *E,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB);

//-- zrot ----------------------------------------------------------------------
void
LAPACK_DECL(zrot)(const INTEGER            *N,
                  DOUBLE_COMPLEX           *CX,
                  const INTEGER            *INCX,
                  DOUBLE_COMPLEX           *CY,
                  const INTEGER            *INCY,
                  const DOUBLE             *C,
                  const DOUBLE_COMPLEX     *S);

//-- zspcon --------------------------------------------------------------------
void
LAPACK_DECL(zspcon)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    const INTEGER            *IPIV,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zspmv ---------------------------------------------------------------------
void
LAPACK_DECL(zspmv)(const char               *UPLO,
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
LAPACK_DECL(zspr)(const char               *UPLO,
                  const INTEGER            *N,
                  const DOUBLE_COMPLEX     *ALPHA,
                  const DOUBLE_COMPLEX     *X,
                  const INTEGER            *INCX,
                  DOUBLE_COMPLEX           *AP);

//-- zsprfs --------------------------------------------------------------------
void
LAPACK_DECL(zsprfs)(const char               *UPLO,
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
LAPACK_DECL(zspsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *AP,
                   INTEGER              *IPIV,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO);

//-- zspsvx --------------------------------------------------------------------
void
LAPACK_DECL(zspsvx)(const char               *FACT,
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
LAPACK_DECL(zsptrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zsptri --------------------------------------------------------------------
void
LAPACK_DECL(zsptri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    const INTEGER    *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zsptrs --------------------------------------------------------------------
void
LAPACK_DECL(zsptrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    const INTEGER            *IPIV,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO);

//-- zstedc --------------------------------------------------------------------
void
LAPACK_DECL(zstedc)(const char       *COMPZ,
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
LAPACK_DECL(zstegr)(const char       *JOBZ,
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
LAPACK_DECL(zstein)(const INTEGER    *N,
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
LAPACK_DECL(zstemr)(const char       *JOBZ,
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
LAPACK_DECL(zsteqr)(const char       *COMPZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO);

//-- zsycon --------------------------------------------------------------------
void
LAPACK_DECL(zsycon)(const char               *UPLO,
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
LAPACK_DECL(zsyconv)(const char       *UPLO,
                     const char       *WAY,
                     const INTEGER    *N,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE_COMPLEX   *WORK,
                     INTEGER          *INFO);

//-- zsyequb -------------------------------------------------------------------
void
LAPACK_DECL(zsyequb)(const char               *UPLO,
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
LAPACK_DECL(zsymv)(const char               *UPLO,
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
LAPACK_DECL(zsyr)(const char               *UPLO,
                  const INTEGER            *N,
                  const DOUBLE_COMPLEX     *ALPHA,
                  const DOUBLE_COMPLEX     *X,
                  const INTEGER            *INCX,
                  DOUBLE_COMPLEX           *A,
                  const INTEGER            *LDA);

//-- zsyrfs --------------------------------------------------------------------
void
LAPACK_DECL(zsyrfs)(const char               *UPLO,
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
LAPACK_DECL(zsyrfsx)(const char               *UPLO,
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
LAPACK_DECL(zsysv)(const char           *UPLO,
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
LAPACK_DECL(zsysvx)(const char               *FACT,
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
LAPACK_DECL(zsysvxx)(const char       *FACT,
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
LAPACK_DECL(zsyswapr)(const char       *UPLO,
                      const INTEGER    *N,
                      DOUBLE_COMPLEX   *A,
                      const INTEGER    *LDA,
                      const INTEGER    *I1,
                      const INTEGER    *I2);

//-- zsytf2 --------------------------------------------------------------------
void
LAPACK_DECL(zsytf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- zsytrf --------------------------------------------------------------------
void
LAPACK_DECL(zsytrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zsytri --------------------------------------------------------------------
void
LAPACK_DECL(zsytri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO);

//-- zsytri2 -------------------------------------------------------------------
void
LAPACK_DECL(zsytri2)(const char       *UPLO,
                     const INTEGER    *N,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE_COMPLEX   *WORK,
                     const INTEGER    *LWORK,
                     INTEGER          *INFO);

//-- zsytri2x ------------------------------------------------------------------
void
LAPACK_DECL(zsytri2x)(const char       *UPLO,
                      const INTEGER    *N,
                      DOUBLE_COMPLEX   *A,
                      const INTEGER    *LDA,
                      const INTEGER    *IPIV,
                      DOUBLE_COMPLEX   *WORK,
                      const INTEGER    *NB,
                      INTEGER          *INFO);

//-- zsytrs --------------------------------------------------------------------
void
LAPACK_DECL(zsytrs)(const char               *UPLO,
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
LAPACK_DECL(zsytrs2)(const char       *UPLO,
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
LAPACK_DECL(ztbcon)(const char               *NORM,
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
LAPACK_DECL(ztbrfs)(const char               *UPLO,
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
LAPACK_DECL(ztbtrs)(const char               *UPLO,
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
LAPACK_DECL(ztfsm)(const char               *TRANSR,
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
LAPACK_DECL(ztftri)(const char       *TRANSR,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    INTEGER          *INFO);

//-- ztfttp --------------------------------------------------------------------
void
LAPACK_DECL(ztfttp)(const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *ARF,
                    DOUBLE_COMPLEX           *AP,
                    INTEGER                  *INFO);

//-- ztfttr --------------------------------------------------------------------
void
LAPACK_DECL(ztfttr)(const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *ARF,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    INTEGER                  *INFO);

//-- ztgevc --------------------------------------------------------------------
void
LAPACK_DECL(ztgevc)(const char               *SIDE,
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
LAPACK_DECL(ztgex2)(const LOGICAL    *WANTQ,
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
LAPACK_DECL(ztgexc)(const LOGICAL    *WANTQ,
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
LAPACK_DECL(ztgsen)(const INTEGER    *IJOB,
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
LAPACK_DECL(ztgsja)(const char       *JOBU,
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
LAPACK_DECL(ztgsna)(const char               *JOB,
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
LAPACK_DECL(ztgsy2)(const char               *TRANS,
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
LAPACK_DECL(ztgsyl)(const char               *TRANS,
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
LAPACK_DECL(ztpcon)(const char               *NORM,
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
LAPACK_DECL(ztprfs)(const char               *UPLO,
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
LAPACK_DECL(ztptri)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    INTEGER          *INFO);

//-- ztptrs --------------------------------------------------------------------
void
LAPACK_DECL(ztptrs)(const char               *UPLO,
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
LAPACK_DECL(ztpttf)(const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE_COMPLEX           *ARF,
                    INTEGER                  *INFO);

//-- ztpttr --------------------------------------------------------------------
void
LAPACK_DECL(ztpttr)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    INTEGER                  *INFO);

//-- ztrcon --------------------------------------------------------------------
void
LAPACK_DECL(ztrcon)(const char               *NORM,
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
LAPACK_DECL(ztrevc)(const char       *SIDE,
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
LAPACK_DECL(ztrexc)(const char       *COMPQ,
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
LAPACK_DECL(ztrrfs)(const char               *UPLO,
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
LAPACK_DECL(ztrsen)(const char       *JOB,
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
LAPACK_DECL(ztrsna)(const char               *JOB,
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
LAPACK_DECL(ztrsyl)(const char               *TRANA,
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
LAPACK_DECL(ztrti2)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- ztrtri --------------------------------------------------------------------
void
LAPACK_DECL(ztrtri)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- ztrtrs --------------------------------------------------------------------
void
LAPACK_DECL(ztrtrs)(const char               *UPLO,
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
LAPACK_DECL(ztrttf)(const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *ARF,
                    INTEGER                  *INFO);

//-- ztrttp --------------------------------------------------------------------
void
LAPACK_DECL(ztrttp)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *AP,
                    INTEGER                  *INFO);

//-- ztzrqf --------------------------------------------------------------------
void
LAPACK_DECL(ztzrqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    INTEGER          *INFO);

//-- ztzrzf --------------------------------------------------------------------
void
LAPACK_DECL(ztzrzf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO);

//-- zunbdb --------------------------------------------------------------------
void
LAPACK_DECL(zunbdb)(const char       *TRANS,
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
LAPACK_DECL(zuncsd)(const char               *JOBU1,
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
LAPACK_DECL(zung2l)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zung2r --------------------------------------------------------------------
void
LAPACK_DECL(zung2r)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zungbr --------------------------------------------------------------------
void
LAPACK_DECL(zungbr)(const char               *VECT,
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
LAPACK_DECL(zunghr)(const INTEGER            *N,
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
LAPACK_DECL(zungl2)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zunglq --------------------------------------------------------------------
void
LAPACK_DECL(zunglq)(const INTEGER            *M,
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
LAPACK_DECL(zungql)(const INTEGER            *M,
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
LAPACK_DECL(zungqr)(const INTEGER            *M,
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
LAPACK_DECL(zungr2)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zungrq --------------------------------------------------------------------
void
LAPACK_DECL(zungrq)(const INTEGER            *M,
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
LAPACK_DECL(zungtr)(const char               *UPLO,
                    const INTEGER            *N,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO);

//-- zunm2l --------------------------------------------------------------------
void
LAPACK_DECL(zunm2l)(const char               *SIDE,
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
LAPACK_DECL(zunm2r)(const char               *SIDE,
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
LAPACK_DECL(zunmbr)(const char               *VECT,
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
LAPACK_DECL(zunmhr)(const char               *SIDE,
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
LAPACK_DECL(zunml2)(const char               *SIDE,
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
LAPACK_DECL(zunmlq)(const char               *SIDE,
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
LAPACK_DECL(zunmql)(const char               *SIDE,
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
LAPACK_DECL(zunmqr)(const char               *SIDE,
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
LAPACK_DECL(zunmr2)(const char               *SIDE,
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
LAPACK_DECL(zunmr3)(const char               *SIDE,
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
LAPACK_DECL(zunmrq)(const char               *SIDE,
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
LAPACK_DECL(zunmrz)(const char               *SIDE,
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
LAPACK_DECL(zunmtr)(const char               *SIDE,
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
LAPACK_DECL(zupgtr)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *Q,
                    const INTEGER            *LDQ,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO);

//-- zupmtr --------------------------------------------------------------------
void
LAPACK_DECL(zupmtr)(const char               *SIDE,
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
LAPACK_DECL(cgetrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    FLOAT_COMPLEX    *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO);

//-- cgetrs --------------------------------------------------------------------
void
LAPACK_DECL(cgetrs)(const char       *TRANS,
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
LAPACK_DECL(clag2z)(const INTEGER    *M,
                    const INTEGER    *N,
                    FLOAT_COMPLEX    *SA,
                    const INTEGER    *LDSA,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- cpotrf --------------------------------------------------------------------
void
LAPACK_DECL(cpotrf)(const char       *UPLO,
                    const INTEGER    *N,
                    FLOAT_COMPLEX    *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- cpotrs --------------------------------------------------------------------
void
LAPACK_DECL(cpotrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    FLOAT_COMPLEX    *A,
                    const INTEGER    *LDA,
                    FLOAT_COMPLEX    *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO);

//-- ilaprec -------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilaprec)(const char   *PREC);

//-- chla_transtype ------------------------------------------------------------
char
LAPACK_DECL(chla_transtype)(const INTEGER    *TRANS);

//-- claswp --------------------------------------------------------------------
void
LAPACK_DECL(claswp)(const INTEGER    *N,
                    FLOAT_COMPLEX    *A,
                    const INTEGER    *LDA,
                    const INTEGER    *K1,
                    const INTEGER    *K2,
                    const INTEGER    *IPIV,
                    const INTEGER    *INCX);

//-- izmax1 --------------------------------------------------------------------
INTEGER
LAPACK_DECL(izmax1)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *CX,
                    const INTEGER            *INCX);

//-- ilazlc --------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilazlc)(const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA);

//-- ilazlr --------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilazlr)(const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA);

//-- cpotf2 --------------------------------------------------------------------
void
LAPACK_DECL(cpotf2)(const char       *UPLO,
                    const INTEGER    *N,
                    FLOAT_COMPLEX    *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO);

//-- clacgv --------------------------------------------------------------------
void
LAPACK_DECL(clacgv)(const INTEGER    *N,
                    FLOAT_COMPLEX    *X,
                    const INTEGER    *INCX);

