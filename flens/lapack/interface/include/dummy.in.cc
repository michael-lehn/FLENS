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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dbbcsd");
    LAPACK_IMPL(dbbcsd)(JOBU1,
                        JOBU2,
                        JOBV1T,
                        JOBV2T,
                        TRANS,
                        M,
                        P,
                        Q,
                        THETA,
                        PHI,
                        U1,
                        LDU1,
                        U2,
                        LDU2,
                        V1T,
                        LDV1T,
                        V2T,
                        LDV2T,
                        B11D,
                        B11E,
                        B12D,
                        B12E,
                        B21D,
                        B21E,
                        B22D,
                        B22E,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dbdsdc");
    LAPACK_IMPL(dbdsdc)(UPLO,
                        COMPQ,
                        N,
                        D,
                        E,
                        U,
                        LDU,
                        VT,
                        LDVT,
                        Q,
                        IQ,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dbdsqr");
    LAPACK_IMPL(dbdsqr)(UPLO,
                        N,
                        NCVT,
                        NRU,
                        NCC,
                        D,
                        E,
                        VT,
                        LDVT,
                        U,
                        LDU,
                        C,
                        LDC,
                        WORK,
                        INFO);
}

//-- ddisna --------------------------------------------------------------------
void
LAPACK_DECL(ddisna)(const char       *JOB,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *D,
                    DOUBLE           *SEP,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("ddisna");
    LAPACK_IMPL(ddisna)(JOB,
                        M,
                        N,
                        D,
                        SEP,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgbbrd");
    LAPACK_IMPL(dgbbrd)(VECT,
                        M,
                        N,
                        NCC,
                        KL,
                        KU,
                        AB,
                        LDAB,
                        D,
                        E,
                        Q,
                        LDQ,
                        PT,
                        LDPT,
                        C,
                        LDC,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgbcon");
    LAPACK_IMPL(dgbcon)(NORM,
                        N,
                        KL,
                        KU,
                        AB,
                        LDAB,
                        IPIV,
                        ANORM,
                        RCOND,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgbequ");
    LAPACK_IMPL(dgbequ)(M,
                        N,
                        KL,
                        KU,
                        AB,
                        LDAB,
                        R,
                        C,
                        ROWCND,
                        COLCND,
                        AMAX,
                        INFO);
}

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
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgbequb");
    LAPACK_IMPL(dgbequb)(M,
                         N,
                         KL,
                         KU,
                         AB,
                         LDAB,
                         R,
                         C,
                         ROWCND,
                         COLCND,
                         AMAX,
                         INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgbrfs");
    LAPACK_IMPL(dgbrfs)(TRANS,
                        N,
                        KL,
                        KU,
                        NRHS,
                        AB,
                        LDAB,
                        AFB,
                        LDAFB,
                        IPIV,
                        B,
                        LDB,
                        X,
                        LDX,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dgbsv");
    LAPACK_IMPL(dgbsv)(N,
                       KL,
                       KU,
                       NRHS,
                       AB,
                       LDAB,
                       IPIV,
                       B,
                       LDB,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgbsvx");
    LAPACK_IMPL(dgbsvx)(FACT,
                        TRANS,
                        N,
                        KL,
                        KU,
                        NRHS,
                        AB,
                        LDAB,
                        AFB,
                        LDAFB,
                        IPIV,
                        EQUED,
                        R,
                        C,
                        B,
                        LDB,
                        X,
                        LDX,
                        RCOND,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dgbtf2 --------------------------------------------------------------------
void
LAPACK_DECL(dgbtf2)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgbtf2");
    LAPACK_IMPL(dgbtf2)(M,
                        N,
                        KL,
                        KU,
                        AB,
                        LDAB,
                        IPIV,
                        INFO);
}

//-- dgbtrf --------------------------------------------------------------------
void
LAPACK_DECL(dgbtrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgbtrf");
    LAPACK_IMPL(dgbtrf)(M,
                        N,
                        KL,
                        KU,
                        AB,
                        LDAB,
                        IPIV,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgbtrs");
    LAPACK_IMPL(dgbtrs)(TRANS,
                        N,
                        KL,
                        KU,
                        NRHS,
                        AB,
                        LDAB,
                        IPIV,
                        B,
                        LDB,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgebak");
    LAPACK_IMPL(dgebak)(JOB,
                        SIDE,
                        N,
                        ILO,
                        IHI,
                        SCALE,
                        M,
                        V,
                        LDV,
                        INFO);
}

//-- dgebal --------------------------------------------------------------------
void
LAPACK_DECL(dgebal)(const char       *JOB,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *ILO,
                    INTEGER          *IHI,
                    DOUBLE           *SCALE,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgebal");
    LAPACK_IMPL(dgebal)(JOB,
                        N,
                        A,
                        LDA,
                        ILO,
                        IHI,
                        SCALE,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgebd2");
    LAPACK_IMPL(dgebd2)(M,
                        N,
                        A,
                        LDA,
                        D,
                        E,
                        TAUQ,
                        TAUP,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgebrd");
    LAPACK_IMPL(dgebrd)(M,
                        N,
                        A,
                        LDA,
                        D,
                        E,
                        TAUQ,
                        TAUP,
                        WORK,
                        LWORK,
                        INFO);
}

//-- dgecon --------------------------------------------------------------------
/*
void
LAPACK_DECL(dgecon)(const char       *NORM,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgecon");
    LAPACK_IMPL(dgecon)(NORM,
                        N,
                        A,
                        LDA,
                        ANORM,
                        RCOND,
                        WORK,
                        IWORK,
                        INFO);
}
*/
//-- dgeequ --------------------------------------------------------------------
/*
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgeequ");
    LAPACK_IMPL(dgeequ)(M,
                        N,
                        A,
                        LDA,
                        R,
                        C,
                        ROWCND,
                        COLCND,
                        AMAX,
                        INFO);
}
*/
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
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgeequb");
    LAPACK_IMPL(dgeequb)(M,
                         N,
                         A,
                         LDA,
                         R,
                         C,
                         ROWCND,
                         COLCND,
                         AMAX,
                         INFO);
}

//-- dgees ---------------------------------------------------------------------
/*
void
LAPACK_DECL(dgees)(const char           *JOBVS,
                   const char           *SORT,
                   const LOGICAL        *SELECT,
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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dgees");
    LAPACK_IMPL(dgees)(JOBVS,
                       SORT,
                       SELECT,
                       N,
                       A,
                       LDA,
                       SDIM,
                       WR,
                       WI,
                       VS,
                       LDVS,
                       WORK,
                       LWORK,
                       BWORK,
                       INFO);
}
*/
//-- dgeesx --------------------------------------------------------------------
/*
void
LAPACK_DECL(dgeesx)(const char       *JOBVS,
                    const char       *SORT,
                    const LOGICAL    *SELECT,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgeesx");
    LAPACK_IMPL(dgeesx)(JOBVS,
                        SORT,
                        SELECT,
                        SENSE,
                        N,
                        A,
                        LDA,
                        SDIM,
                        WR,
                        WI,
                        VS,
                        LDVS,
                        RCONDE,
                        RCONDV,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        BWORK,
                        INFO);
}
*/
//-- dgeev ---------------------------------------------------------------------
/*
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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dgeev");
    LAPACK_IMPL(dgeev)(JOBVL,
                       JOBVR,
                       N,
                       A,
                       LDA,
                       WR,
                       WI,
                       VL,
                       LDVL,
                       VR,
                       LDVR,
                       WORK,
                       LWORK,
                       INFO);
}
*/
//-- dgeevx --------------------------------------------------------------------
/*
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgeevx");
    LAPACK_IMPL(dgeevx)(BALANC,
                        JOBVL,
                        JOBVR,
                        SENSE,
                        N,
                        A,
                        LDA,
                        WR,
                        WI,
                        VL,
                        LDVL,
                        VR,
                        LDVR,
                        ILO,
                        IHI,
                        SCALE,
                        ABNRM,
                        RCONDE,
                        RCONDV,
                        WORK,
                        LWORK,
                        IWORK,
                        INFO);
}
*/
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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dgegs");
    LAPACK_IMPL(dgegs)(JOBVSL,
                       JOBVSR,
                       N,
                       A,
                       LDA,
                       B,
                       LDB,
                       ALPHAR,
                       ALPHAI,
                       BETA,
                       VSL,
                       LDVSL,
                       VSR,
                       LDVSR,
                       WORK,
                       LWORK,
                       INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dgegv");
    LAPACK_IMPL(dgegv)(JOBVL,
                       JOBVR,
                       N,
                       A,
                       LDA,
                       B,
                       LDB,
                       ALPHAR,
                       ALPHAI,
                       BETA,
                       VL,
                       LDVL,
                       VR,
                       LDVR,
                       WORK,
                       LWORK,
                       INFO);
}

//-- dgehd2 --------------------------------------------------------------------
void
LAPACK_DECL(dgehd2)(const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgehd2");
    LAPACK_IMPL(dgehd2)(N,
                        ILO,
                        IHI,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgehrd");
    LAPACK_IMPL(dgehrd)(N,
                        ILO,
                        IHI,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

//-- dgejsv --------------------------------------------------------------------
/*
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgejsv");
    LAPACK_IMPL(dgejsv)(JOBA,
                        JOBU,
                        JOBV,
                        JOBR,
                        JOBT,
                        JOBP,
                        M,
                        N,
                        A,
                        LDA,
                        SVA,
                        U,
                        LDU,
                        V,
                        LDV,
                        WORK,
                        LWORK,
                        IWORK,
                        INFO);
}
*/
//-- dgelq2 --------------------------------------------------------------------
void
LAPACK_DECL(dgelq2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgelq2");
    LAPACK_IMPL(dgelq2)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

//-- dgelqf --------------------------------------------------------------------
/*
void
LAPACK_DECL(dgelqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgelqf");
    LAPACK_IMPL(dgelqf)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}
*/
//-- dgels ---------------------------------------------------------------------
/*
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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dgels");
    LAPACK_IMPL(dgels)(TRANS,
                       M,
                       N,
                       NRHS,
                       A,
                       LDA,
                       B,
                       LDB,
                       WORK,
                       LWORK,
                       INFO);
}
*/
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgelsd");
    LAPACK_IMPL(dgelsd)(M,
                        N,
                        NRHS,
                        A,
                        LDA,
                        B,
                        LDB,
                        S,
                        RCOND,
                        RANK,
                        WORK,
                        LWORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgelss");
    LAPACK_IMPL(dgelss)(M,
                        N,
                        NRHS,
                        A,
                        LDA,
                        B,
                        LDB,
                        S,
                        RCOND,
                        RANK,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgelsx");
    LAPACK_IMPL(dgelsx)(M,
                        N,
                        NRHS,
                        A,
                        LDA,
                        B,
                        LDB,
                        JPVT,
                        RCOND,
                        RANK,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgelsy");
    LAPACK_IMPL(dgelsy)(M,
                        N,
                        NRHS,
                        A,
                        LDA,
                        B,
                        LDB,
                        JPVT,
                        RCOND,
                        RANK,
                        WORK,
                        LWORK,
                        INFO);
}

//-- dgeql2 --------------------------------------------------------------------
void
LAPACK_DECL(dgeql2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgeql2");
    LAPACK_IMPL(dgeql2)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

//-- dgeqlf --------------------------------------------------------------------
void
LAPACK_DECL(dgeqlf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgeqlf");
    LAPACK_IMPL(dgeqlf)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

//-- dgeqp3 --------------------------------------------------------------------
/*
void
LAPACK_DECL(dgeqp3)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *JPVT,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgeqp3");
    LAPACK_IMPL(dgeqp3)(M,
                        N,
                        A,
                        LDA,
                        JPVT,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}
*/

//-- dgeqpf --------------------------------------------------------------------
void
LAPACK_DECL(dgeqpf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *JPVT,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgeqpf");
    LAPACK_IMPL(dgeqpf)(M,
                        N,
                        A,
                        LDA,
                        JPVT,
                        TAU,
                        WORK,
                        INFO);
}

//-- dgeqr2 --------------------------------------------------------------------
void
LAPACK_DECL(dgeqr2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgeqr2");
    LAPACK_IMPL(dgeqr2)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

//-- dgeqr2p -------------------------------------------------------------------
void
LAPACK_DECL(dgeqr2p)(const INTEGER    *M,
                     const INTEGER    *N,
                     DOUBLE           *A,
                     const INTEGER    *LDA,
                     DOUBLE           *TAU,
                     DOUBLE           *WORK,
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgeqr2p");
    LAPACK_IMPL(dgeqr2p)(M,
                         N,
                         A,
                         LDA,
                         TAU,
                         WORK,
                         INFO);
}

//-- dgeqrf --------------------------------------------------------------------
/*
void
LAPACK_DECL(dgeqrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgeqrf");
    LAPACK_IMPL(dgeqrf)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}
*/
//-- dgeqrfp -------------------------------------------------------------------
void
LAPACK_DECL(dgeqrfp)(const INTEGER    *M,
                     const INTEGER    *N,
                     DOUBLE           *A,
                     const INTEGER    *LDA,
                     DOUBLE           *TAU,
                     DOUBLE           *WORK,
                     const INTEGER    *LWORK,
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgeqrfp");
    LAPACK_IMPL(dgeqrfp)(M,
                         N,
                         A,
                         LDA,
                         TAU,
                         WORK,
                         LWORK,
                         INFO);
}

//-- dgerfs --------------------------------------------------------------------
/*
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgerfs");
    LAPACK_IMPL(dgerfs)(TRANS,
                        N,
                        NRHS,
                        A,
                        LDA,
                        AF,
                        LDAF,
                        IPIV,
                        B,
                        LDB,
                        X,
                        LDX,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}
*/
//-- dgerq2 --------------------------------------------------------------------
void
LAPACK_DECL(dgerq2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgerq2");
    LAPACK_IMPL(dgerq2)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

//-- dgerqf --------------------------------------------------------------------
void
LAPACK_DECL(dgerqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgerqf");
    LAPACK_IMPL(dgerqf)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

//-- dgesc2 --------------------------------------------------------------------
void
LAPACK_DECL(dgesc2)(const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *RHS,
                    const INTEGER    *IPIV,
                    const INTEGER    *JPIV,
                    DOUBLE           *SCALE)
{
    DEBUG_LAPACK_STUB("dgesc2");
    LAPACK_IMPL(dgesc2)(N,
                        A,
                        LDA,
                        RHS,
                        IPIV,
                        JPIV,
                        SCALE);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgesdd");
    LAPACK_IMPL(dgesdd)(JOBZ,
                        M,
                        N,
                        A,
                        LDA,
                        S,
                        U,
                        LDU,
                        VT,
                        LDVT,
                        WORK,
                        LWORK,
                        IWORK,
                        INFO);
}

//-- dgesv ---------------------------------------------------------------------
/*
void
LAPACK_DECL(dgesv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   INTEGER              *IPIV,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dgesv");
    LAPACK_IMPL(dgesv)(N,
                       NRHS,
                       A,
                       LDA,
                       IPIV,
                       B,
                       LDB,
                       INFO);
}
*/
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgesvd");
    LAPACK_IMPL(dgesvd)(JOBU,
                        JOBVT,
                        M,
                        N,
                        A,
                        LDA,
                        S,
                        U,
                        LDU,
                        VT,
                        LDVT,
                        WORK,
                        LWORK,
                        INFO);
}

//-- dgesvj --------------------------------------------------------------------
/*
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgesvj");
    LAPACK_IMPL(dgesvj)(JOBA,
                        JOBU,
                        JOBV,
                        M,
                        N,
                        A,
                        LDA,
                        SVA,
                        MV,
                        V,
                        LDV,
                        WORK,
                        LWORK,
                        INFO);
}
*/
//-- dgesvx --------------------------------------------------------------------
/*
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgesvx");
    LAPACK_IMPL(dgesvx)(FACT,
                        TRANS,
                        N,
                        NRHS,
                        A,
                        LDA,
                        AF,
                        LDAF,
                        IPIV,
                        EQUED,
                        R,
                        C,
                        B,
                        LDB,
                        X,
                        LDX,
                        RCOND,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}
*/
//-- dgetc2 --------------------------------------------------------------------
void
LAPACK_DECL(dgetc2)(const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *JPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgetc2");
    LAPACK_IMPL(dgetc2)(N,
                        A,
                        LDA,
                        IPIV,
                        JPIV,
                        INFO);
}

//-- dgetf2 --------------------------------------------------------------------
void
LAPACK_DECL(dgetf2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgetf2");
    LAPACK_IMPL(dgetf2)(M,
                        N,
                        A,
                        LDA,
                        IPIV,
                        INFO);
}

//-- dgetrf --------------------------------------------------------------------
/*
void
LAPACK_DECL(dgetrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgetrf");
    LAPACK_IMPL(dgetrf)(M,
                        N,
                        A,
                        LDA,
                        IPIV,
                        INFO);
}
*/
//-- dgetri --------------------------------------------------------------------
/*
void
LAPACK_DECL(dgetri)(const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgetri");
    LAPACK_IMPL(dgetri)(N,
                        A,
                        LDA,
                        IPIV,
                        WORK,
                        LWORK,
                        INFO);
}
*/
//-- dgetrs --------------------------------------------------------------------
/*
void
LAPACK_DECL(dgetrs)(const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgetrs");
    LAPACK_IMPL(dgetrs)(TRANS,
                        N,
                        NRHS,
                        A,
                        LDA,
                        IPIV,
                        B,
                        LDB,
                        INFO);
}
*/
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dggbak");
    LAPACK_IMPL(dggbak)(JOB,
                        SIDE,
                        N,
                        ILO,
                        IHI,
                        LSCALE,
                        RSCALE,
                        M,
                        V,
                        LDV,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dggbal");
    LAPACK_IMPL(dggbal)(JOB,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        ILO,
                        IHI,
                        LSCALE,
                        RSCALE,
                        WORK,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dgges");
    LAPACK_IMPL(dgges)(JOBVSL,
                       JOBVSR,
                       SORT,
                       SELCTG,
                       N,
                       A,
                       LDA,
                       B,
                       LDB,
                       SDIM,
                       ALPHAR,
                       ALPHAI,
                       BETA,
                       VSL,
                       LDVSL,
                       VSR,
                       LDVSR,
                       WORK,
                       LWORK,
                       BWORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dggesx");
    LAPACK_IMPL(dggesx)(JOBVSL,
                        JOBVSR,
                        SORT,
                        SELCTG,
                        SENSE,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        SDIM,
                        ALPHAR,
                        ALPHAI,
                        BETA,
                        VSL,
                        LDVSL,
                        VSR,
                        LDVSR,
                        RCONDE,
                        RCONDV,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        BWORK,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dggev");
    LAPACK_IMPL(dggev)(JOBVL,
                       JOBVR,
                       N,
                       A,
                       LDA,
                       B,
                       LDB,
                       ALPHAR,
                       ALPHAI,
                       BETA,
                       VL,
                       LDVL,
                       VR,
                       LDVR,
                       WORK,
                       LWORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dggevx");
    LAPACK_IMPL(dggevx)(BALANC,
                        JOBVL,
                        JOBVR,
                        SENSE,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        ALPHAR,
                        ALPHAI,
                        BETA,
                        VL,
                        LDVL,
                        VR,
                        LDVR,
                        ILO,
                        IHI,
                        LSCALE,
                        RSCALE,
                        ABNRM,
                        BBNRM,
                        RCONDE,
                        RCONDV,
                        WORK,
                        LWORK,
                        IWORK,
                        BWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dggglm");
    LAPACK_IMPL(dggglm)(N,
                        M,
                        P,
                        A,
                        LDA,
                        B,
                        LDB,
                        D,
                        X,
                        Y,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgghrd");
    LAPACK_IMPL(dgghrd)(COMPQ,
                        COMPZ,
                        N,
                        ILO,
                        IHI,
                        A,
                        LDA,
                        B,
                        LDB,
                        Q,
                        LDQ,
                        Z,
                        LDZ,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgglse");
    LAPACK_IMPL(dgglse)(M,
                        N,
                        P,
                        A,
                        LDA,
                        B,
                        LDB,
                        C,
                        D,
                        X,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dggqrf");
    LAPACK_IMPL(dggqrf)(N,
                        M,
                        P,
                        A,
                        LDA,
                        TAUA,
                        B,
                        LDB,
                        TAUB,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dggrqf");
    LAPACK_IMPL(dggrqf)(M,
                        P,
                        N,
                        A,
                        LDA,
                        TAUA,
                        B,
                        LDB,
                        TAUB,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dggsvd");
    LAPACK_IMPL(dggsvd)(JOBU,
                        JOBV,
                        JOBQ,
                        M,
                        N,
                        P,
                        K,
                        L,
                        A,
                        LDA,
                        B,
                        LDB,
                        ALPHA,
                        BETA,
                        U,
                        LDU,
                        V,
                        LDV,
                        Q,
                        LDQ,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dggsvp");
    LAPACK_IMPL(dggsvp)(JOBU,
                        JOBV,
                        JOBQ,
                        M,
                        P,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        TOLA,
                        TOLB,
                        K,
                        L,
                        U,
                        LDU,
                        V,
                        LDV,
                        Q,
                        LDQ,
                        IWORK,
                        TAU,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgsvj0");
    LAPACK_IMPL(dgsvj0)(JOBV,
                        M,
                        N,
                        A,
                        LDA,
                        D,
                        SVA,
                        MV,
                        V,
                        LDV,
                        EPS,
                        SFMIN,
                        TOL,
                        NSWEEP,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgsvj1");
    LAPACK_IMPL(dgsvj1)(JOBV,
                        M,
                        N,
                        N1,
                        A,
                        LDA,
                        D,
                        SVA,
                        MV,
                        V,
                        LDV,
                        EPS,
                        SFMIN,
                        TOL,
                        NSWEEP,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgtcon");
    LAPACK_IMPL(dgtcon)(NORM,
                        N,
                        DL,
                        D,
                        DU,
                        DU2,
                        IPIV,
                        ANORM,
                        RCOND,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgtrfs");
    LAPACK_IMPL(dgtrfs)(TRANS,
                        N,
                        NRHS,
                        DL,
                        D,
                        DU,
                        DLF,
                        DF,
                        DUF,
                        DU2,
                        IPIV,
                        B,
                        LDB,
                        X,
                        LDX,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dgtsv ---------------------------------------------------------------------
void
LAPACK_DECL(dgtsv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *DL,
                   DOUBLE               *D,
                   DOUBLE               *DU,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dgtsv");
    LAPACK_IMPL(dgtsv)(N,
                       NRHS,
                       DL,
                       D,
                       DU,
                       B,
                       LDB,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgtsvx");
    LAPACK_IMPL(dgtsvx)(FACT,
                        TRANS,
                        N,
                        NRHS,
                        DL,
                        D,
                        DU,
                        DLF,
                        DF,
                        DUF,
                        DU2,
                        IPIV,
                        B,
                        LDB,
                        X,
                        LDX,
                        RCOND,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dgttrf --------------------------------------------------------------------
void
LAPACK_DECL(dgttrf)(const INTEGER    *N,
                    DOUBLE           *DL,
                    DOUBLE           *D,
                    DOUBLE           *DU,
                    DOUBLE           *DU2,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgttrf");
    LAPACK_IMPL(dgttrf)(N,
                        DL,
                        D,
                        DU,
                        DU2,
                        IPIV,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dgttrs");
    LAPACK_IMPL(dgttrs)(TRANS,
                        N,
                        NRHS,
                        DL,
                        D,
                        DU,
                        DU2,
                        IPIV,
                        B,
                        LDB,
                        INFO);
}

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
                    const INTEGER    *LDB)
{
    DEBUG_LAPACK_STUB("dgtts2");
    LAPACK_IMPL(dgtts2)(ITRANS,
                        N,
                        NRHS,
                        DL,
                        D,
                        DU,
                        DU2,
                        IPIV,
                        B,
                        LDB);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dhgeqz");
    LAPACK_IMPL(dhgeqz)(JOB,
                        COMPQ,
                        COMPZ,
                        N,
                        ILO,
                        IHI,
                        H,
                        LDH,
                        T,
                        LDT,
                        ALPHAR,
                        ALPHAI,
                        BETA,
                        Q,
                        LDQ,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dhsein");
    LAPACK_IMPL(dhsein)(SIDE,
                        EIGSRC,
                        INITV,
                        SELECT,
                        N,
                        H,
                        LDH,
                        WR,
                        WI,
                        VL,
                        LDVL,
                        VR,
                        LDVR,
                        MM,
                        M,
                        WORK,
                        IFAILL,
                        IFAILR,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dhseqr");
    LAPACK_IMPL(dhseqr)(JOB,
                        COMPZ,
                        N,
                        ILO,
                        IHI,
                        H,
                        LDH,
                        WR,
                        WI,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        INFO);
}

//-- disnan --------------------------------------------------------------------
LOGICAL
LAPACK_DECL(disnan)(const DOUBLE     *DIN)
{
    DEBUG_LAPACK_STUB("disnan");
    return LAPACK_IMPL(disnan)(DIN);
}

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
                       const INTEGER        *INCY)
{
    DEBUG_LAPACK_STUB("dla_gbamv");
    LAPACK_IMPL(dla_gbamv)(TRANS,
                           M,
                           N,
                           KL,
                           KU,
                           ALPHA,
                           AB,
                           LDAB,
                           X,
                           INCX,
                           BETA,
                           Y,
                           INCY);
}

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
                         const INTEGER    *IWORK)
{
    DEBUG_LAPACK_STUB("dla_gbrcond");
    return LAPACK_IMPL(dla_gbrcond)(TRANS,
                                    N,
                                    KL,
                                    KU,
                                    AB,
                                    LDAB,
                                    AFB,
                                    LDAFB,
                                    IPIV,
                                    CMODE,
                                    C,
                                    INFO,
                                    WORK,
                                    IWORK);
}

//-- dla_gbrpvgrw --------------------------------------------------------------
DOUBLE
LAPACK_DECL(dla_gbrpvgrw)(const INTEGER    *N,
                          const INTEGER    *KL,
                          const INTEGER    *KU,
                          const INTEGER    *NCOLS,
                          const DOUBLE     *AB,
                          const INTEGER    *LDAB,
                          const DOUBLE     *AFB,
                          const INTEGER    *LDAFB)
{
    DEBUG_LAPACK_STUB("dla_gbrpvgrw");
    return LAPACK_IMPL(dla_gbrpvgrw)(N,
                                     KL,
                                     KU,
                                     NCOLS,
                                     AB,
                                     LDAB,
                                     AFB,
                                     LDAFB);
}

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
                       const INTEGER        *INCY)
{
    DEBUG_LAPACK_STUB("dla_geamv");
    LAPACK_IMPL(dla_geamv)(TRANS,
                           M,
                           N,
                           ALPHA,
                           A,
                           LDA,
                           X,
                           INCX,
                           BETA,
                           Y,
                           INCY);
}

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
                         const INTEGER    *IWORK)
{
    DEBUG_LAPACK_STUB("dla_gercond");
    return LAPACK_IMPL(dla_gercond)(TRANS,
                                    N,
                                    A,
                                    LDA,
                                    AF,
                                    LDAF,
                                    IPIV,
                                    CMODE,
                                    C,
                                    INFO,
                                    WORK,
                                    IWORK);
}

//-- dla_lin_berr --------------------------------------------------------------
void
LAPACK_DECL(dla_lin_berr)(const INTEGER    *N,
                          const INTEGER    *NZ,
                          const INTEGER    *NRHS,
                          const DOUBLE     *RES,
                          const DOUBLE     *AYB,
                          DOUBLE           *BERR)
{
    DEBUG_LAPACK_STUB("dla_lin_berr");
    LAPACK_IMPL(dla_lin_berr)(N,
                              NZ,
                              NRHS,
                              RES,
                              AYB,
                              BERR);
}

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
                         const INTEGER    *IWORK)
{
    DEBUG_LAPACK_STUB("dla_porcond");
    return LAPACK_IMPL(dla_porcond)(UPLO,
                                    N,
                                    A,
                                    LDA,
                                    AF,
                                    LDAF,
                                    CMODE,
                                    C,
                                    INFO,
                                    WORK,
                                    IWORK);
}

//-- dla_porpvgrw --------------------------------------------------------------
DOUBLE
LAPACK_DECL(dla_porpvgrw)(const char       *UPLO,
                          const INTEGER    *NCOLS,
                          const DOUBLE     *A,
                          const INTEGER    *LDA,
                          const DOUBLE     *AF,
                          const INTEGER    *LDAF,
                          const DOUBLE     *WORK)
{
    DEBUG_LAPACK_STUB("dla_porpvgrw");
    return LAPACK_IMPL(dla_porpvgrw)(UPLO,
                                     NCOLS,
                                     A,
                                     LDA,
                                     AF,
                                     LDAF,
                                     WORK);
}

//-- dla_rpvgrw ----------------------------------------------------------------
DOUBLE
LAPACK_DECL(dla_rpvgrw)(const INTEGER    *N,
                        const INTEGER    *NCOLS,
                        const DOUBLE     *A,
                        const INTEGER    *LDA,
                        const DOUBLE     *AF,
                        const INTEGER    *LDAF)
{
    DEBUG_LAPACK_STUB("dla_rpvgrw");
    return LAPACK_IMPL(dla_rpvgrw)(N,
                                   NCOLS,
                                   A,
                                   LDA,
                                   AF,
                                   LDAF);
}

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
                       const INTEGER        *INCY)
{
    DEBUG_LAPACK_STUB("dla_syamv");
    LAPACK_IMPL(dla_syamv)(UPLO,
                           N,
                           ALPHA,
                           A,
                           LDA,
                           X,
                           INCX,
                           BETA,
                           Y,
                           INCY);
}

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
                         const INTEGER    *IWORK)
{
    DEBUG_LAPACK_STUB("dla_syrcond");
    return LAPACK_IMPL(dla_syrcond)(UPLO,
                                    N,
                                    A,
                                    LDA,
                                    AF,
                                    LDAF,
                                    IPIV,
                                    CMODE,
                                    C,
                                    INFO,
                                    WORK,
                                    IWORK);
}

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
                          const DOUBLE     *WORK)
{
    DEBUG_LAPACK_STUB("dla_syrpvgrw");
    return LAPACK_IMPL(dla_syrpvgrw)(UPLO,
                                     N,
                                     INFO,
                                     A,
                                     LDA,
                                     AF,
                                     LDAF,
                                     IPIV,
                                     WORK);
}

//-- dla_wwaddw ----------------------------------------------------------------
void
LAPACK_DECL(dla_wwaddw)(const INTEGER    *N,
                        DOUBLE           *X,
                        DOUBLE           *Y,
                        const DOUBLE     *W)
{
    DEBUG_LAPACK_STUB("dla_wwaddw");
    LAPACK_IMPL(dla_wwaddw)(N,
                            X,
                            Y,
                            W);
}

//-- dlabad --------------------------------------------------------------------
void
LAPACK_DECL(dlabad)(DOUBLE   *SMALL,
                    DOUBLE   *LARGE)
{
    DEBUG_LAPACK_STUB("dlabad");
    LAPACK_IMPL(dlabad)(SMALL,
                        LARGE);
}

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
                    const INTEGER    *LDY)
{
    DEBUG_LAPACK_STUB("dlabrd");
    LAPACK_IMPL(dlabrd)(M,
                        N,
                        NB,
                        A,
                        LDA,
                        D,
                        E,
                        TAUQ,
                        TAUP,
                        X,
                        LDX,
                        Y,
                        LDY);
}

//-- dlacn2 --------------------------------------------------------------------
void
LAPACK_DECL(dlacn2)(const INTEGER    *N,
                    DOUBLE           *V,
                    DOUBLE           *X,
                    INTEGER          *ISGN,
                    DOUBLE           *EST,
                    INTEGER          *KASE,
                    INTEGER          *ISAVE)
{
    DEBUG_LAPACK_STUB("dlacn2");
    LAPACK_IMPL(dlacn2)(N,
                        V,
                        X,
                        ISGN,
                        EST,
                        KASE,
                        ISAVE);
}

//-- dlacon --------------------------------------------------------------------
void
LAPACK_DECL(dlacon)(const INTEGER    *N,
                    DOUBLE           *V,
                    DOUBLE           *X,
                    INTEGER          *ISGN,
                    DOUBLE           *EST,
                    INTEGER          *KASE)
{
    DEBUG_LAPACK_STUB("dlacon");
    LAPACK_IMPL(dlacon)(N,
                        V,
                        X,
                        ISGN,
                        EST,
                        KASE);
}

//-- dlacpy --------------------------------------------------------------------
void
LAPACK_DECL(dlacpy)(const char       *UPLO,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB)
{
    DEBUG_LAPACK_STUB("dlacpy");
    LAPACK_IMPL(dlacpy)(UPLO,
                        M,
                        N,
                        A,
                        LDA,
                        B,
                        LDB);
}

//-- dladiv --------------------------------------------------------------------
void
LAPACK_DECL(dladiv)(const DOUBLE     *A,
                    const DOUBLE     *B,
                    const DOUBLE     *C,
                    const DOUBLE     *D,
                    DOUBLE           *P,
                    DOUBLE           *Q)
{
    DEBUG_LAPACK_STUB("dladiv");
    LAPACK_IMPL(dladiv)(A,
                        B,
                        C,
                        D,
                        P,
                        Q);
}

//-- dlae2 ---------------------------------------------------------------------
void
LAPACK_DECL(dlae2)(const DOUBLE     *A,
                   const DOUBLE     *B,
                   const DOUBLE     *C,
                   DOUBLE           *RT1,
                   DOUBLE           *RT2)
{
    DEBUG_LAPACK_STUB("dlae2");
    LAPACK_IMPL(dlae2)(A,
                       B,
                       C,
                       RT1,
                       RT2);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaebz");
    LAPACK_IMPL(dlaebz)(IJOB,
                        NITMAX,
                        N,
                        MMAX,
                        MINP,
                        NBMIN,
                        ABSTOL,
                        RELTOL,
                        PIVMIN,
                        D,
                        E,
                        E2,
                        NVAL,
                        AB,
                        C,
                        MOUT,
                        NAB,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaed0");
    LAPACK_IMPL(dlaed0)(ICOMPQ,
                        QSIZ,
                        N,
                        D,
                        E,
                        Q,
                        LDQ,
                        QSTORE,
                        LDQS,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaed1");
    LAPACK_IMPL(dlaed1)(N,
                        D,
                        Q,
                        LDQ,
                        INDXQ,
                        RHO,
                        CUTPNT,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaed2");
    LAPACK_IMPL(dlaed2)(K,
                        N,
                        N1,
                        D,
                        Q,
                        LDQ,
                        INDXQ,
                        RHO,
                        Z,
                        DLAMDA,
                        W,
                        Q2,
                        INDX,
                        INDXC,
                        INDXP,
                        COLTYP,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaed3");
    LAPACK_IMPL(dlaed3)(K,
                        N,
                        N1,
                        D,
                        Q,
                        LDQ,
                        RHO,
                        DLAMDA,
                        Q2,
                        INDX,
                        CTOT,
                        W,
                        S,
                        INFO);
}

//-- dlaed4 --------------------------------------------------------------------
void
LAPACK_DECL(dlaed4)(const INTEGER    *N,
                    const INTEGER    *I,
                    const DOUBLE     *D,
                    const DOUBLE     *Z,
                    DOUBLE           *DELTA,
                    const DOUBLE     *RHO,
                    DOUBLE           *DLAM,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaed4");
    LAPACK_IMPL(dlaed4)(N,
                        I,
                        D,
                        Z,
                        DELTA,
                        RHO,
                        DLAM,
                        INFO);
}

//-- dlaed5 --------------------------------------------------------------------
void
LAPACK_DECL(dlaed5)(const INTEGER    *I,
                    const DOUBLE     *D,
                    const DOUBLE     *Z,
                    DOUBLE           *DELTA,
                    const DOUBLE     *RHO,
                    DOUBLE           *DLAM)
{
    DEBUG_LAPACK_STUB("dlaed5");
    LAPACK_IMPL(dlaed5)(I,
                        D,
                        Z,
                        DELTA,
                        RHO,
                        DLAM);
}

//-- dlaed6 --------------------------------------------------------------------
void
LAPACK_DECL(dlaed6)(const INTEGER    *KNITER,
                    const LOGICAL    *ORGATI,
                    const DOUBLE     *RHO,
                    const DOUBLE     *D,
                    const DOUBLE     *Z,
                    const DOUBLE     *FINIT,
                    DOUBLE           *TAU,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaed6");
    LAPACK_IMPL(dlaed6)(KNITER,
                        ORGATI,
                        RHO,
                        D,
                        Z,
                        FINIT,
                        TAU,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaed7");
    LAPACK_IMPL(dlaed7)(ICOMPQ,
                        N,
                        QSIZ,
                        TLVLS,
                        CURLVL,
                        CURPBM,
                        D,
                        Q,
                        LDQ,
                        INDXQ,
                        RHO,
                        CUTPNT,
                        QSTORE,
                        QPTR,
                        PRMPTR,
                        PERM,
                        GIVPTR,
                        GIVCOL,
                        GIVNUM,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaed8");
    LAPACK_IMPL(dlaed8)(ICOMPQ,
                        K,
                        N,
                        QSIZ,
                        D,
                        Q,
                        LDQ,
                        INDXQ,
                        RHO,
                        CUTPNT,
                        Z,
                        DLAMDA,
                        Q2,
                        LDQ2,
                        W,
                        PERM,
                        GIVPTR,
                        GIVCOL,
                        GIVNUM,
                        INDXP,
                        INDX,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaed9");
    LAPACK_IMPL(dlaed9)(K,
                        KSTART,
                        KSTOP,
                        N,
                        D,
                        Q,
                        LDQ,
                        RHO,
                        DLAMDA,
                        W,
                        S,
                        LDS,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaeda");
    LAPACK_IMPL(dlaeda)(N,
                        TLVLS,
                        CURLVL,
                        CURPBM,
                        PRMPTR,
                        PERM,
                        GIVPTR,
                        GIVCOL,
                        GIVNUM,
                        Q,
                        QPTR,
                        Z,
                        ZTEMP,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaein");
    LAPACK_IMPL(dlaein)(RIGHTV,
                        NOINIT,
                        N,
                        H,
                        LDH,
                        WR,
                        WI,
                        VR,
                        VI,
                        B,
                        LDB,
                        WORK,
                        EPS3,
                        SMLNUM,
                        BIGNUM,
                        INFO);
}

//-- dlaev2 --------------------------------------------------------------------
void
LAPACK_DECL(dlaev2)(const DOUBLE     *A,
                    const DOUBLE     *B,
                    const DOUBLE     *C,
                    DOUBLE           *RT1,
                    DOUBLE           *RT2,
                    DOUBLE           *CS1,
                    DOUBLE           *SN1)
{
    DEBUG_LAPACK_STUB("dlaev2");
    LAPACK_IMPL(dlaev2)(A,
                        B,
                        C,
                        RT1,
                        RT2,
                        CS1,
                        SN1);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaexc");
    LAPACK_IMPL(dlaexc)(WANTQ,
                        N,
                        T,
                        LDT,
                        Q,
                        LDQ,
                        J1,
                        N1,
                        N2,
                        WORK,
                        INFO);
}

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
                   DOUBLE               *WI)
{
    DEBUG_LAPACK_STUB("dlag2");
    LAPACK_IMPL(dlag2)(A,
                       LDA,
                       B,
                       LDB,
                       SAFMIN,
                       SCALE1,
                       SCALE2,
                       WR1,
                       WR2,
                       WI);
}

//-- dlag2s --------------------------------------------------------------------
void
LAPACK_DECL(dlag2s)(const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    FLOAT            *SA,
                    const INTEGER    *LDSA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlag2s");
    LAPACK_IMPL(dlag2s)(M,
                        N,
                        A,
                        LDA,
                        SA,
                        LDSA,
                        INFO);
}

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
                    DOUBLE           *SNQ)
{
    DEBUG_LAPACK_STUB("dlags2");
    LAPACK_IMPL(dlags2)(UPPER,
                        A1,
                        A2,
                        A3,
                        B1,
                        B2,
                        B3,
                        CSU,
                        SNU,
                        CSV,
                        SNV,
                        CSQ,
                        SNQ);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlagtf");
    LAPACK_IMPL(dlagtf)(N,
                        A,
                        LAMBDA,
                        B,
                        C,
                        TOL,
                        D,
                        IN,
                        INFO);
}

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
                    const INTEGER    *LDB)
{
    DEBUG_LAPACK_STUB("dlagtm");
    LAPACK_IMPL(dlagtm)(TRANS,
                        N,
                        NRHS,
                        ALPHA,
                        DL,
                        D,
                        DU,
                        X,
                        LDX,
                        BETA,
                        B,
                        LDB);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlagts");
    LAPACK_IMPL(dlagts)(JOB,
                        N,
                        A,
                        B,
                        C,
                        D,
                        IN,
                        Y,
                        TOL,
                        INFO);
}

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
                    DOUBLE           *SNR)
{
    DEBUG_LAPACK_STUB("dlagv2");
    LAPACK_IMPL(dlagv2)(A,
                        LDA,
                        B,
                        LDB,
                        ALPHAR,
                        ALPHAI,
                        BETA,
                        CSL,
                        SNL,
                        CSR,
                        SNR);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlahqr");
    LAPACK_IMPL(dlahqr)(WANTT,
                        WANTZ,
                        N,
                        ILO,
                        IHI,
                        H,
                        LDH,
                        WR,
                        WI,
                        ILOZ,
                        IHIZ,
                        Z,
                        LDZ,
                        INFO);
}

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
                    const INTEGER    *LDY)
{
    DEBUG_LAPACK_STUB("dlahr2");
    LAPACK_IMPL(dlahr2)(N,
                        K,
                        NB,
                        A,
                        LDA,
                        TAU,
                        T,
                        LDT,
                        Y,
                        LDY);
}

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
                    const INTEGER    *LDY)
{
    DEBUG_LAPACK_STUB("dlahrd");
    LAPACK_IMPL(dlahrd)(N,
                        K,
                        NB,
                        A,
                        LDA,
                        TAU,
                        T,
                        LDT,
                        Y,
                        LDY);
}

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
                    DOUBLE           *C)
{
    DEBUG_LAPACK_STUB("dlaic1");
    LAPACK_IMPL(dlaic1)(JOB,
                        J,
                        X,
                        SEST,
                        W,
                        GAMMA,
                        SESTPR,
                        S,
                        C);
}

//-- dlaisnan ------------------------------------------------------------------
LOGICAL
LAPACK_DECL(dlaisnan)(const DOUBLE     *DIN1,
                      const DOUBLE     *DIN2)
{
    DEBUG_LAPACK_STUB("dlaisnan");
    return LAPACK_IMPL(dlaisnan)(DIN1,
                                 DIN2);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaln2");
    LAPACK_IMPL(dlaln2)(LTRANS,
                        NA,
                        NW,
                        SMIN,
                        CA,
                        A,
                        LDA,
                        D1,
                        D2,
                        B,
                        LDB,
                        WR,
                        WI,
                        X,
                        LDX,
                        SCALE,
                        XNORM,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlals0");
    LAPACK_IMPL(dlals0)(ICOMPQ,
                        NL,
                        NR,
                        SQRE,
                        NRHS,
                        B,
                        LDB,
                        BX,
                        LDBX,
                        PERM,
                        GIVPTR,
                        GIVCOL,
                        LDGCOL,
                        GIVNUM,
                        LDGNUM,
                        POLES,
                        DIFL,
                        DIFR,
                        Z,
                        K,
                        C,
                        S,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlalsa");
    LAPACK_IMPL(dlalsa)(ICOMPQ,
                        SMLSIZ,
                        N,
                        NRHS,
                        B,
                        LDB,
                        BX,
                        LDBX,
                        U,
                        LDU,
                        VT,
                        K,
                        DIFL,
                        DIFR,
                        Z,
                        POLES,
                        GIVPTR,
                        GIVCOL,
                        LDGCOL,
                        PERM,
                        GIVNUM,
                        C,
                        S,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlalsd");
    LAPACK_IMPL(dlalsd)(UPLO,
                        SMLSIZ,
                        N,
                        NRHS,
                        D,
                        E,
                        B,
                        LDB,
                        RCOND,
                        RANK,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dlamch --------------------------------------------------------------------
/*
DOUBLE
LAPACK_DECL(dlamch)(const char   *CMACH)
{
    DEBUG_LAPACK_STUB("dlamch");
    return LAPACK_IMPL(dlamch)(CMACH);
}
*/
//-- dlamrg --------------------------------------------------------------------
void
LAPACK_DECL(dlamrg)(const INTEGER    *N1,
                    const INTEGER    *N2,
                    const DOUBLE     *A,
                    const INTEGER    *DTRD1,
                    const INTEGER    *DTRD2,
                    INTEGER          *INDEX)
{
    DEBUG_LAPACK_STUB("dlamrg");
    LAPACK_IMPL(dlamrg)(N1,
                        N2,
                        A,
                        DTRD1,
                        DTRD2,
                        INDEX);
}

//-- dlaneg --------------------------------------------------------------------
INTEGER
LAPACK_DECL(dlaneg)(const INTEGER    *N,
                    const DOUBLE     *D,
                    const DOUBLE     *LLD,
                    const DOUBLE     *SIGMA,
                    const DOUBLE     *PIVMIN,
                    const INTEGER    *R)
{
    DEBUG_LAPACK_STUB("dlaneg");
    return LAPACK_IMPL(dlaneg)(N,
                               D,
                               LLD,
                               SIGMA,
                               PIVMIN,
                               R);
}

//-- dlangb --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlangb)(const char       *NORM,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlangb");
    return LAPACK_IMPL(dlangb)(NORM,
                               N,
                               KL,
                               KU,
                               AB,
                               LDAB,
                               WORK);
}

//-- dlange --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlange)(const char       *NORM,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlange");
    return LAPACK_IMPL(dlange)(NORM,
                               M,
                               N,
                               A,
                               LDA,
                               WORK);
}

//-- dlangt --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlangt)(const char       *NORM,
                    const INTEGER    *N,
                    const DOUBLE     *DL,
                    const DOUBLE     *D,
                    const DOUBLE     *DU)
{
    DEBUG_LAPACK_STUB("dlangt");
    return LAPACK_IMPL(dlangt)(NORM,
                               N,
                               DL,
                               D,
                               DU);
}

//-- dlanhs --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlanhs)(const char       *NORM,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlanhs");
    return LAPACK_IMPL(dlanhs)(NORM,
                               N,
                               A,
                               LDA,
                               WORK);
}

//-- dlansb --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlansb)(const char       *NORM,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlansb");
    return LAPACK_IMPL(dlansb)(NORM,
                               UPLO,
                               N,
                               K,
                               AB,
                               LDAB,
                               WORK);
}

//-- dlansf --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlansf)(const char       *NORM,
                    const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlansf");
    return LAPACK_IMPL(dlansf)(NORM,
                               TRANSR,
                               UPLO,
                               N,
                               A,
                               WORK);
}

//-- dlansp --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlansp)(const char       *NORM,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlansp");
    return LAPACK_IMPL(dlansp)(NORM,
                               UPLO,
                               N,
                               AP,
                               WORK);
}

//-- dlanst --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlanst)(const char       *NORM,
                    const INTEGER    *N,
                    const DOUBLE     *D,
                    const DOUBLE     *E)
{
    DEBUG_LAPACK_STUB("dlanst");
    return LAPACK_IMPL(dlanst)(NORM,
                               N,
                               D,
                               E);
}

//-- dlansy --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlansy)(const char       *NORM,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlansy");
    return LAPACK_IMPL(dlansy)(NORM,
                               UPLO,
                               N,
                               A,
                               LDA,
                               WORK);
}

//-- dlantb --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlantb)(const char       *NORM,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    const DOUBLE     *AB,
                    const INTEGER    *LDAB,
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlantb");
    return LAPACK_IMPL(dlantb)(NORM,
                               UPLO,
                               DIAG,
                               N,
                               K,
                               AB,
                               LDAB,
                               WORK);
}

//-- dlantp --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlantp)(const char       *NORM,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlantp");
    return LAPACK_IMPL(dlantp)(NORM,
                               UPLO,
                               DIAG,
                               N,
                               AP,
                               WORK);
}

//-- dlantr --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlantr)(const char       *NORM,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlantr");
    return LAPACK_IMPL(dlantr)(NORM,
                               UPLO,
                               DIAG,
                               M,
                               N,
                               A,
                               LDA,
                               WORK);
}

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
                    DOUBLE   *SN)
{
    DEBUG_LAPACK_STUB("dlanv2");
    LAPACK_IMPL(dlanv2)(A,
                        B,
                        C,
                        D,
                        RT1R,
                        RT1I,
                        RT2R,
                        RT2I,
                        CS,
                        SN);
}

//-- dlapll --------------------------------------------------------------------
void
LAPACK_DECL(dlapll)(const INTEGER    *N,
                    DOUBLE           *X,
                    const INTEGER    *INCX,
                    DOUBLE           *Y,
                    const INTEGER    *INCY,
                    DOUBLE           *SSMIN)
{
    DEBUG_LAPACK_STUB("dlapll");
    LAPACK_IMPL(dlapll)(N,
                        X,
                        INCX,
                        Y,
                        INCY,
                        SSMIN);
}

//-- dlapmr --------------------------------------------------------------------
void
LAPACK_DECL(dlapmr)(const LOGICAL    *FORWRD,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    INTEGER          *K)
{
    DEBUG_LAPACK_STUB("dlapmr");
    LAPACK_IMPL(dlapmr)(FORWRD,
                        M,
                        N,
                        X,
                        LDX,
                        K);
}

//-- dlapmt --------------------------------------------------------------------
void
LAPACK_DECL(dlapmt)(const LOGICAL    *FORWRD,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    INTEGER          *K)
{
    DEBUG_LAPACK_STUB("dlapmt");
    LAPACK_IMPL(dlapmt)(FORWRD,
                        M,
                        N,
                        X,
                        LDX,
                        K);
}

//-- dlapy2 --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlapy2)(const DOUBLE     *X,
                    const DOUBLE     *Y)
{
    DEBUG_LAPACK_STUB("dlapy2");
    return LAPACK_IMPL(dlapy2)(X,
                               Y);
}

//-- dlapy3 --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dlapy3)(const DOUBLE     *X,
                    const DOUBLE     *Y,
                    const DOUBLE     *Z)
{
    DEBUG_LAPACK_STUB("dlapy3");
    return LAPACK_IMPL(dlapy3)(X,
                               Y,
                               Z);
}

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
                    char             *EQUED)
{
    DEBUG_LAPACK_STUB("dlaqgb");
    LAPACK_IMPL(dlaqgb)(M,
                        N,
                        KL,
                        KU,
                        AB,
                        LDAB,
                        R,
                        C,
                        ROWCND,
                        COLCND,
                        AMAX,
                        EQUED);
}

//-- dlaqge --------------------------------------------------------------------
/*
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
                    char             *EQUED)
{
    DEBUG_LAPACK_STUB("dlaqge");
    LAPACK_IMPL(dlaqge)(M,
                        N,
                        A,
                        LDA,
                        R,
                        C,
                        ROWCND,
                        COLCND,
                        AMAX,
                        EQUED);
}
*/
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
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlaqp2");
    LAPACK_IMPL(dlaqp2)(M,
                        N,
                        OFFSET,
                        A,
                        LDA,
                        JPVT,
                        TAU,
                        VN1,
                        VN2,
                        WORK);
}

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
                    const INTEGER    *LDF)
{
    DEBUG_LAPACK_STUB("dlaqps");
    LAPACK_IMPL(dlaqps)(M,
                        N,
                        OFFSET,
                        NB,
                        KB,
                        A,
                        LDA,
                        JPVT,
                        TAU,
                        VN1,
                        VN2,
                        AUXV,
                        F,
                        LDF);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaqr0");
    LAPACK_IMPL(dlaqr0)(WANTT,
                        WANTZ,
                        N,
                        ILO,
                        IHI,
                        H,
                        LDH,
                        WR,
                        WI,
                        ILOZ,
                        IHIZ,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        INFO);
}

//-- dlaqr1 --------------------------------------------------------------------
void
LAPACK_DECL(dlaqr1)(const INTEGER    *N,
                    const DOUBLE     *H,
                    const INTEGER    *LDH,
                    const DOUBLE     *SR1,
                    const DOUBLE     *SI1,
                    const DOUBLE     *SR2,
                    const DOUBLE     *SI2,
                    DOUBLE           *V)
{
    DEBUG_LAPACK_STUB("dlaqr1");
    LAPACK_IMPL(dlaqr1)(N,
                        H,
                        LDH,
                        SR1,
                        SI1,
                        SR2,
                        SI2,
                        V);
}

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
                    const INTEGER    *LWORK)
{
    DEBUG_LAPACK_STUB("dlaqr2");
    LAPACK_IMPL(dlaqr2)(WANTT,
                        WANTZ,
                        N,
                        KTOP,
                        KBOT,
                        NW,
                        H,
                        LDH,
                        ILOZ,
                        IHIZ,
                        Z,
                        LDZ,
                        NS,
                        ND,
                        SR,
                        SI,
                        V,
                        LDV,
                        NH,
                        T,
                        LDT,
                        NV,
                        WV,
                        LDWV,
                        WORK,
                        LWORK);
}

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
                    const INTEGER    *LWORK)
{
    DEBUG_LAPACK_STUB("dlaqr3");
    LAPACK_IMPL(dlaqr3)(WANTT,
                        WANTZ,
                        N,
                        KTOP,
                        KBOT,
                        NW,
                        H,
                        LDH,
                        ILOZ,
                        IHIZ,
                        Z,
                        LDZ,
                        NS,
                        ND,
                        SR,
                        SI,
                        V,
                        LDV,
                        NH,
                        T,
                        LDT,
                        NV,
                        WV,
                        LDWV,
                        WORK,
                        LWORK);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaqr4");
    LAPACK_IMPL(dlaqr4)(WANTT,
                        WANTZ,
                        N,
                        ILO,
                        IHI,
                        H,
                        LDH,
                        WR,
                        WI,
                        ILOZ,
                        IHIZ,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    const INTEGER    *LDWH)
{
    DEBUG_LAPACK_STUB("dlaqr5");
    LAPACK_IMPL(dlaqr5)(WANTT,
                        WANTZ,
                        KACC22,
                        N,
                        KTOP,
                        KBOT,
                        NSHFTS,
                        SR,
                        SI,
                        H,
                        LDH,
                        ILOZ,
                        IHIZ,
                        Z,
                        LDZ,
                        V,
                        LDV,
                        U,
                        LDU,
                        NV,
                        WV,
                        LDWV,
                        NH,
                        WH,
                        LDWH);
}

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
                    char             *EQUED)
{
    DEBUG_LAPACK_STUB("dlaqsb");
    LAPACK_IMPL(dlaqsb)(UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        S,
                        SCOND,
                        AMAX,
                        EQUED);
}

//-- dlaqsp --------------------------------------------------------------------
void
LAPACK_DECL(dlaqsp)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED)
{
    DEBUG_LAPACK_STUB("dlaqsp");
    LAPACK_IMPL(dlaqsp)(UPLO,
                        N,
                        AP,
                        S,
                        SCOND,
                        AMAX,
                        EQUED);
}

//-- dlaqsy --------------------------------------------------------------------
void
LAPACK_DECL(dlaqsy)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED)
{
    DEBUG_LAPACK_STUB("dlaqsy");
    LAPACK_IMPL(dlaqsy)(UPLO,
                        N,
                        A,
                        LDA,
                        S,
                        SCOND,
                        AMAX,
                        EQUED);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlaqtr");
    LAPACK_IMPL(dlaqtr)(LTRAN,
                        LREAL,
                        N,
                        T,
                        LDT,
                        B,
                        W,
                        SCALE,
                        X,
                        WORK,
                        INFO);
}

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
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlar1v");
    LAPACK_IMPL(dlar1v)(N,
                        B1,
                        BN,
                        LAMBDA,
                        D,
                        L,
                        LD,
                        LLD,
                        PIVMIN,
                        GAPTOL,
                        Z,
                        WANTNC,
                        NEGCNT,
                        ZTZ,
                        MINGMA,
                        R,
                        ISUPPZ,
                        NRMINV,
                        RESID,
                        RQCORR,
                        WORK);
}

//-- dlar2v --------------------------------------------------------------------
void
LAPACK_DECL(dlar2v)(const INTEGER    *N,
                    DOUBLE           *X,
                    DOUBLE           *Y,
                    DOUBLE           *Z,
                    const INTEGER    *INCX,
                    const DOUBLE     *C,
                    const DOUBLE     *S,
                    const INTEGER    *INCC)
{
    DEBUG_LAPACK_STUB("dlar2v");
    LAPACK_IMPL(dlar2v)(N,
                        X,
                        Y,
                        Z,
                        INCX,
                        C,
                        S,
                        INCC);
}

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
                   DOUBLE               *WORK)
{
    DEBUG_LAPACK_STUB("dlarf");
    LAPACK_IMPL(dlarf)(SIDE,
                       M,
                       N,
                       V,
                       INCV,
                       TAU,
                       C,
                       LDC,
                       WORK);
}

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
                    const INTEGER    *LDWORK)
{
    DEBUG_LAPACK_STUB("dlarfb");
    LAPACK_IMPL(dlarfb)(SIDE,
                        TRANS,
                        DIRECT,
                        STOREV,
                        M,
                        N,
                        K,
                        V,
                        LDV,
                        T,
                        LDT,
                        C,
                        LDC,
                        WORK,
                        LDWORK);
}

//-- dlarfg --------------------------------------------------------------------
void
LAPACK_DECL(dlarfg)(const INTEGER    *N,
                    DOUBLE           *ALPHA,
                    DOUBLE           *X,
                    const INTEGER    *INCX,
                    DOUBLE           *TAU)
{
    DEBUG_LAPACK_STUB("dlarfg");
    LAPACK_IMPL(dlarfg)(N,
                        ALPHA,
                        X,
                        INCX,
                        TAU);
}

//-- dlarfgp -------------------------------------------------------------------
void
LAPACK_DECL(dlarfgp)(const INTEGER    *N,
                     DOUBLE           *ALPHA,
                     DOUBLE           *X,
                     const INTEGER    *INCX,
                     DOUBLE           *TAU)
{
    DEBUG_LAPACK_STUB("dlarfgp");
    LAPACK_IMPL(dlarfgp)(N,
                         ALPHA,
                         X,
                         INCX,
                         TAU);
}

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
                    const INTEGER    *LDT)
{
    DEBUG_LAPACK_STUB("dlarft");
    LAPACK_IMPL(dlarft)(DIRECT,
                        STOREV,
                        N,
                        K,
                        V,
                        LDV,
                        TAU,
                        T,
                        LDT);
}

//-- dlarfx --------------------------------------------------------------------
void
LAPACK_DECL(dlarfx)(const char       *SIDE,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *V,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlarfx");
    LAPACK_IMPL(dlarfx)(SIDE,
                        M,
                        N,
                        V,
                        TAU,
                        C,
                        LDC,
                        WORK);
}

//-- dlargv --------------------------------------------------------------------
void
LAPACK_DECL(dlargv)(const INTEGER    *N,
                    DOUBLE           *X,
                    const INTEGER    *INCX,
                    DOUBLE           *Y,
                    const INTEGER    *INCY,
                    DOUBLE           *C,
                    const INTEGER    *INCC)
{
    DEBUG_LAPACK_STUB("dlargv");
    LAPACK_IMPL(dlargv)(N,
                        X,
                        INCX,
                        Y,
                        INCY,
                        C,
                        INCC);
}

//-- dlarnv --------------------------------------------------------------------
void
LAPACK_DECL(dlarnv)(const INTEGER    *IDIST,
                    INTEGER          *ISEED,
                    const INTEGER    *N,
                    DOUBLE           *X)
{
    DEBUG_LAPACK_STUB("dlarnv");
    LAPACK_IMPL(dlarnv)(IDIST,
                        ISEED,
                        N,
                        X);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlarra");
    LAPACK_IMPL(dlarra)(N,
                        D,
                        E,
                        E2,
                        SPLTOL,
                        TNRM,
                        NSPLIT,
                        ISPLIT,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlarrb");
    LAPACK_IMPL(dlarrb)(N,
                        D,
                        LLD,
                        IFIRST,
                        ILAST,
                        RTOL1,
                        RTOL2,
                        OFFSET,
                        W,
                        WGAP,
                        WERR,
                        WORK,
                        IWORK,
                        PIVMIN,
                        SPDIAM,
                        TWIST,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlarrc");
    LAPACK_IMPL(dlarrc)(JOBT,
                        N,
                        VL,
                        VU,
                        D,
                        E,
                        PIVMIN,
                        EIGCNT,
                        LCNT,
                        RCNT,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlarrd");
    LAPACK_IMPL(dlarrd)(RANGE,
                        ORDER,
                        N,
                        VL,
                        VU,
                        IL,
                        IU,
                        GERS,
                        RELTOL,
                        D,
                        E,
                        E2,
                        PIVMIN,
                        NSPLIT,
                        ISPLIT,
                        M,
                        W,
                        WERR,
                        WL,
                        WU,
                        IBLOCK,
                        INDEXW,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlarre");
    LAPACK_IMPL(dlarre)(RANGE,
                        N,
                        VL,
                        VU,
                        IL,
                        IU,
                        D,
                        E,
                        E2,
                        RTOL1,
                        RTOL2,
                        SPLTOL,
                        NSPLIT,
                        ISPLIT,
                        M,
                        W,
                        WERR,
                        WGAP,
                        IBLOCK,
                        INDEXW,
                        GERS,
                        PIVMIN,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlarrf");
    LAPACK_IMPL(dlarrf)(N,
                        D,
                        L,
                        LD,
                        CLSTRT,
                        CLEND,
                        W,
                        WGAP,
                        WERR,
                        SPDIAM,
                        CLGAPL,
                        CLGAPR,
                        PIVMIN,
                        SIGMA,
                        DPLUS,
                        LPLUS,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlarrj");
    LAPACK_IMPL(dlarrj)(N,
                        D,
                        E2,
                        IFIRST,
                        ILAST,
                        RTOL,
                        OFFSET,
                        W,
                        WERR,
                        WORK,
                        IWORK,
                        PIVMIN,
                        SPDIAM,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlarrk");
    LAPACK_IMPL(dlarrk)(N,
                        IW,
                        GL,
                        GU,
                        D,
                        E2,
                        PIVMIN,
                        RELTOL,
                        W,
                        WERR,
                        INFO);
}

//-- dlarrr --------------------------------------------------------------------
void
LAPACK_DECL(dlarrr)(const INTEGER    *N,
                    const DOUBLE     *D,
                    DOUBLE           *E,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlarrr");
    LAPACK_IMPL(dlarrr)(N,
                        D,
                        E,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlarrv");
    LAPACK_IMPL(dlarrv)(N,
                        VL,
                        VU,
                        D,
                        L,
                        PIVMIN,
                        ISPLIT,
                        M,
                        DOL,
                        DOU,
                        MINRGP,
                        RTOL1,
                        RTOL2,
                        W,
                        WERR,
                        WGAP,
                        IBLOCK,
                        INDEXW,
                        GERS,
                        Z,
                        LDZ,
                        ISUPPZ,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dlarscl2 ------------------------------------------------------------------
void
LAPACK_DECL(dlarscl2)(const INTEGER    *M,
                      const INTEGER    *N,
                      const DOUBLE     *D,
                      DOUBLE           *X,
                      const INTEGER    *LDX)
{
    DEBUG_LAPACK_STUB("dlarscl2");
    LAPACK_IMPL(dlarscl2)(M,
                          N,
                          D,
                          X,
                          LDX);
}

//-- dlartg --------------------------------------------------------------------
void
LAPACK_DECL(dlartg)(const DOUBLE     *F,
                    const DOUBLE     *G,
                    DOUBLE           *CS,
                    DOUBLE           *SN,
                    DOUBLE           *R)
{
    DEBUG_LAPACK_STUB("dlartg");
    LAPACK_IMPL(dlartg)(F,
                        G,
                        CS,
                        SN,
                        R);
}

//-- dlartgp -------------------------------------------------------------------
void
LAPACK_DECL(dlartgp)(const DOUBLE     *F,
                     const DOUBLE     *G,
                     DOUBLE           *CS,
                     DOUBLE           *SN,
                     DOUBLE           *R)
{
    DEBUG_LAPACK_STUB("dlartgp");
    LAPACK_IMPL(dlartgp)(F,
                         G,
                         CS,
                         SN,
                         R);
}

//-- dlartgs -------------------------------------------------------------------
void
LAPACK_DECL(dlartgs)(const DOUBLE     *X,
                     const DOUBLE     *Y,
                     const DOUBLE     *SIGMA,
                     DOUBLE           *CS,
                     DOUBLE           *SN)
{
    DEBUG_LAPACK_STUB("dlartgs");
    LAPACK_IMPL(dlartgs)(X,
                         Y,
                         SIGMA,
                         CS,
                         SN);
}

//-- dlartv --------------------------------------------------------------------
void
LAPACK_DECL(dlartv)(const INTEGER    *N,
                    DOUBLE           *X,
                    const INTEGER    *INCX,
                    DOUBLE           *Y,
                    const INTEGER    *INCY,
                    const DOUBLE     *C,
                    const DOUBLE     *S,
                    const INTEGER    *INCC)
{
    DEBUG_LAPACK_STUB("dlartv");
    LAPACK_IMPL(dlartv)(N,
                        X,
                        INCX,
                        Y,
                        INCY,
                        C,
                        S,
                        INCC);
}

//-- dlaruv --------------------------------------------------------------------
void
LAPACK_DECL(dlaruv)(INTEGER          *ISEED,
                    const INTEGER    *N,
                    DOUBLE           *X)
{
    DEBUG_LAPACK_STUB("dlaruv");
    LAPACK_IMPL(dlaruv)(ISEED,
                        N,
                        X);
}

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
                   DOUBLE               *WORK)
{
    DEBUG_LAPACK_STUB("dlarz");
    LAPACK_IMPL(dlarz)(SIDE,
                       M,
                       N,
                       L,
                       V,
                       INCV,
                       TAU,
                       C,
                       LDC,
                       WORK);
}

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
                    const INTEGER    *LDWORK)
{
    DEBUG_LAPACK_STUB("dlarzb");
    LAPACK_IMPL(dlarzb)(SIDE,
                        TRANS,
                        DIRECT,
                        STOREV,
                        M,
                        N,
                        K,
                        L,
                        V,
                        LDV,
                        T,
                        LDT,
                        C,
                        LDC,
                        WORK,
                        LDWORK);
}

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
                    const INTEGER    *LDT)
{
    DEBUG_LAPACK_STUB("dlarzt");
    LAPACK_IMPL(dlarzt)(DIRECT,
                        STOREV,
                        N,
                        K,
                        V,
                        LDV,
                        TAU,
                        T,
                        LDT);
}

//-- dlas2 ---------------------------------------------------------------------
void
LAPACK_DECL(dlas2)(const DOUBLE     *F,
                   const DOUBLE     *G,
                   const DOUBLE     *H,
                   DOUBLE           *SSMIN,
                   DOUBLE           *SSMAX)
{
    DEBUG_LAPACK_STUB("dlas2");
    LAPACK_IMPL(dlas2)(F,
                       G,
                       H,
                       SSMIN,
                       SSMAX);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlascl");
    LAPACK_IMPL(dlascl)(TYPE,
                        KL,
                        KU,
                        CFROM,
                        CTO,
                        M,
                        N,
                        A,
                        LDA,
                        INFO);
}

//-- dlascl2 -------------------------------------------------------------------
void
LAPACK_DECL(dlascl2)(const INTEGER    *M,
                     const INTEGER    *N,
                     const DOUBLE     *D,
                     DOUBLE           *X,
                     const INTEGER    *LDX)
{
    DEBUG_LAPACK_STUB("dlascl2");
    LAPACK_IMPL(dlascl2)(M,
                         N,
                         D,
                         X,
                         LDX);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlasd0");
    LAPACK_IMPL(dlasd0)(N,
                        SQRE,
                        D,
                        E,
                        U,
                        LDU,
                        VT,
                        LDVT,
                        SMLSIZ,
                        IWORK,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlasd1");
    LAPACK_IMPL(dlasd1)(NL,
                        NR,
                        SQRE,
                        D,
                        ALPHA,
                        BETA,
                        U,
                        LDU,
                        VT,
                        LDVT,
                        IDXQ,
                        IWORK,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlasd2");
    LAPACK_IMPL(dlasd2)(NL,
                        NR,
                        SQRE,
                        K,
                        D,
                        Z,
                        ALPHA,
                        BETA,
                        U,
                        LDU,
                        VT,
                        LDVT,
                        DSIGMA,
                        U2,
                        LDU2,
                        VT2,
                        LDVT2,
                        IDXP,
                        IDX,
                        IDXC,
                        IDXQ,
                        COLTYP,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlasd3");
    LAPACK_IMPL(dlasd3)(NL,
                        NR,
                        SQRE,
                        K,
                        D,
                        Q,
                        LDQ,
                        DSIGMA,
                        U,
                        LDU,
                        U2,
                        LDU2,
                        VT,
                        LDVT,
                        VT2,
                        LDVT2,
                        IDXC,
                        CTOT,
                        Z,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlasd4");
    LAPACK_IMPL(dlasd4)(N,
                        I,
                        D,
                        Z,
                        DELTA,
                        RHO,
                        SIGMA,
                        WORK,
                        INFO);
}

//-- dlasd5 --------------------------------------------------------------------
void
LAPACK_DECL(dlasd5)(const INTEGER    *I,
                    const DOUBLE     *D,
                    const DOUBLE     *Z,
                    DOUBLE           *DELTA,
                    const DOUBLE     *RHO,
                    DOUBLE           *DSIGMA,
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlasd5");
    LAPACK_IMPL(dlasd5)(I,
                        D,
                        Z,
                        DELTA,
                        RHO,
                        DSIGMA,
                        WORK);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlasd6");
    LAPACK_IMPL(dlasd6)(ICOMPQ,
                        NL,
                        NR,
                        SQRE,
                        D,
                        VF,
                        VL,
                        ALPHA,
                        BETA,
                        IDXQ,
                        PERM,
                        GIVPTR,
                        GIVCOL,
                        LDGCOL,
                        GIVNUM,
                        LDGNUM,
                        POLES,
                        DIFL,
                        DIFR,
                        Z,
                        K,
                        C,
                        S,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlasd7");
    LAPACK_IMPL(dlasd7)(ICOMPQ,
                        NL,
                        NR,
                        SQRE,
                        K,
                        D,
                        Z,
                        ZW,
                        VF,
                        VFW,
                        VL,
                        VLW,
                        ALPHA,
                        BETA,
                        DSIGMA,
                        IDX,
                        IDXP,
                        IDXQ,
                        PERM,
                        GIVPTR,
                        GIVCOL,
                        LDGCOL,
                        GIVNUM,
                        LDGNUM,
                        C,
                        S,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlasd8");
    LAPACK_IMPL(dlasd8)(ICOMPQ,
                        K,
                        D,
                        Z,
                        VF,
                        VL,
                        DIFL,
                        DIFR,
                        LDDIFR,
                        DSIGMA,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlasda");
    LAPACK_IMPL(dlasda)(ICOMPQ,
                        SMLSIZ,
                        N,
                        SQRE,
                        D,
                        E,
                        U,
                        LDU,
                        VT,
                        K,
                        DIFL,
                        DIFR,
                        Z,
                        POLES,
                        GIVPTR,
                        GIVCOL,
                        LDGCOL,
                        PERM,
                        GIVNUM,
                        C,
                        S,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlasdq");
    LAPACK_IMPL(dlasdq)(UPLO,
                        SQRE,
                        N,
                        NCVT,
                        NRU,
                        NCC,
                        D,
                        E,
                        VT,
                        LDVT,
                        U,
                        LDU,
                        C,
                        LDC,
                        WORK,
                        INFO);
}

//-- dlasdt --------------------------------------------------------------------
void
LAPACK_DECL(dlasdt)(const INTEGER    *N,
                    INTEGER          *LVL,
                    INTEGER          *ND,
                    INTEGER          *INODE,
                    INTEGER          *NDIML,
                    INTEGER          *NDIMR,
                    const INTEGER    *MSUB)
{
    DEBUG_LAPACK_STUB("dlasdt");
    LAPACK_IMPL(dlasdt)(N,
                        LVL,
                        ND,
                        INODE,
                        NDIML,
                        NDIMR,
                        MSUB);
}

//-- dlaset --------------------------------------------------------------------
void
LAPACK_DECL(dlaset)(const char       *UPLO,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *ALPHA,
                    const DOUBLE     *BETA,
                    DOUBLE           *A,
                    const INTEGER    *LDA)
{
    DEBUG_LAPACK_STUB("dlaset");
    LAPACK_IMPL(dlaset)(UPLO,
                        M,
                        N,
                        ALPHA,
                        BETA,
                        A,
                        LDA);
}

//-- dlasq1 --------------------------------------------------------------------
void
LAPACK_DECL(dlasq1)(const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlasq1");
    LAPACK_IMPL(dlasq1)(N,
                        D,
                        E,
                        WORK,
                        INFO);
}

//-- dlasq2 --------------------------------------------------------------------
void
LAPACK_DECL(dlasq2)(const INTEGER    *N,
                    DOUBLE           *Z,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlasq2");
    LAPACK_IMPL(dlasq2)(N,
                        Z,
                        INFO);
}

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
                    DOUBLE           *TAU)
{
    DEBUG_LAPACK_STUB("dlasq3");
    LAPACK_IMPL(dlasq3)(I0,
                        N0,
                        Z,
                        PP,
                        DMIN,
                        SIGMA,
                        DESIG,
                        QMAX,
                        NFAIL,
                        ITER,
                        NDIV,
                        IEEE,
                        TTYPE,
                        DMIN1,
                        DMIN2,
                        DN,
                        DN1,
                        DN2,
                        G,
                        TAU);
}

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
                    DOUBLE           *G)
{
    DEBUG_LAPACK_STUB("dlasq4");
    LAPACK_IMPL(dlasq4)(I0,
                        N0,
                        Z,
                        PP,
                        N0IN,
                        DMIN,
                        DMIN1,
                        DMIN2,
                        DN,
                        DN1,
                        DN2,
                        TAU,
                        TTYPE,
                        G);
}

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
                    const LOGICAL    *IEEE)
{
    DEBUG_LAPACK_STUB("dlasq5");
    LAPACK_IMPL(dlasq5)(I0,
                        N0,
                        Z,
                        PP,
                        TAU,
                        DMIN,
                        DMIN1,
                        DMIN2,
                        DN,
                        DNM1,
                        DNM2,
                        IEEE);
}

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
                    DOUBLE           *DNM2)
{
    DEBUG_LAPACK_STUB("dlasq6");
    LAPACK_IMPL(dlasq6)(I0,
                        N0,
                        Z,
                        PP,
                        DMIN,
                        DMIN1,
                        DMIN2,
                        DN,
                        DNM1,
                        DNM2);
}

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
                   const INTEGER        *LDA)
{
    DEBUG_LAPACK_STUB("dlasr");
    LAPACK_IMPL(dlasr)(SIDE,
                       PIVOT,
                       DIRECT,
                       M,
                       N,
                       C,
                       S,
                       A,
                       LDA);
}

//-- dlasrt --------------------------------------------------------------------
void
LAPACK_DECL(dlasrt)(const char       *ID,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlasrt");
    LAPACK_IMPL(dlasrt)(ID,
                        N,
                        D,
                        INFO);
}

//-- dlassq --------------------------------------------------------------------
void
LAPACK_DECL(dlassq)(const INTEGER    *N,
                    const DOUBLE     *X,
                    const INTEGER    *INCX,
                    DOUBLE           *SCALE,
                    DOUBLE           *SUMSQ)
{
    DEBUG_LAPACK_STUB("dlassq");
    LAPACK_IMPL(dlassq)(N,
                        X,
                        INCX,
                        SCALE,
                        SUMSQ);
}

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
                    DOUBLE           *CSL)
{
    DEBUG_LAPACK_STUB("dlasv2");
    LAPACK_IMPL(dlasv2)(F,
                        G,
                        H,
                        SSMIN,
                        SSMAX,
                        SNR,
                        CSR,
                        SNL,
                        CSL);
}

//-- dlaswp --------------------------------------------------------------------
void
LAPACK_DECL(dlaswp)(const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const INTEGER    *K1,
                    const INTEGER    *K2,
                    const INTEGER    *IPIV,
                    const INTEGER    *INCX)
{
    DEBUG_LAPACK_STUB("dlaswp");
    LAPACK_IMPL(dlaswp)(N,
                        A,
                        LDA,
                        K1,
                        K2,
                        IPIV,
                        INCX);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlasy2");
    LAPACK_IMPL(dlasy2)(LTRANL,
                        LTRANR,
                        ISGN,
                        N1,
                        N2,
                        TL,
                        LDTL,
                        TR,
                        LDTR,
                        B,
                        LDB,
                        SCALE,
                        X,
                        LDX,
                        XNORM,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlasyf");
    LAPACK_IMPL(dlasyf)(UPLO,
                        N,
                        NB,
                        KB,
                        A,
                        LDA,
                        IPIV,
                        W,
                        LDW,
                        INFO);
}

//-- dlat2s --------------------------------------------------------------------
void
LAPACK_DECL(dlat2s)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    FLOAT            *SA,
                    const INTEGER    *LDSA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlat2s");
    LAPACK_IMPL(dlat2s)(UPLO,
                        N,
                        A,
                        LDA,
                        SA,
                        LDSA,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlatbs");
    LAPACK_IMPL(dlatbs)(UPLO,
                        TRANS,
                        DIAG,
                        NORMIN,
                        N,
                        KD,
                        AB,
                        LDAB,
                        X,
                        SCALE,
                        CNORM,
                        INFO);
}

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
                    const INTEGER    *JPIV)
{
    DEBUG_LAPACK_STUB("dlatdf");
    LAPACK_IMPL(dlatdf)(IJOB,
                        N,
                        Z,
                        LDZ,
                        RHS,
                        RDSUM,
                        RDSCAL,
                        IPIV,
                        JPIV);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlatps");
    LAPACK_IMPL(dlatps)(UPLO,
                        TRANS,
                        DIAG,
                        NORMIN,
                        N,
                        AP,
                        X,
                        SCALE,
                        CNORM,
                        INFO);
}

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
                    const INTEGER    *LDW)
{
    DEBUG_LAPACK_STUB("dlatrd");
    LAPACK_IMPL(dlatrd)(UPLO,
                        N,
                        NB,
                        A,
                        LDA,
                        E,
                        TAU,
                        W,
                        LDW);
}

//-- dlatrs --------------------------------------------------------------------
/*
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlatrs");
    LAPACK_IMPL(dlatrs)(UPLO,
                        TRANS,
                        DIAG,
                        NORMIN,
                        N,
                        A,
                        LDA,
                        X,
                        SCALE,
                        CNORM,
                        INFO);
}
*/
//-- dlatrz --------------------------------------------------------------------
void
LAPACK_DECL(dlatrz)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *L,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlatrz");
    LAPACK_IMPL(dlatrz)(M,
                        N,
                        L,
                        A,
                        LDA,
                        TAU,
                        WORK);
}

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
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("dlatzm");
    LAPACK_IMPL(dlatzm)(SIDE,
                        M,
                        N,
                        V,
                        INCV,
                        TAU,
                        C1,
                        C2,
                        LDC,
                        WORK);
}

//-- dlauu2 --------------------------------------------------------------------
void
LAPACK_DECL(dlauu2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlauu2");
    LAPACK_IMPL(dlauu2)(UPLO,
                        N,
                        A,
                        LDA,
                        INFO);
}

//-- dlauum --------------------------------------------------------------------
void
LAPACK_DECL(dlauum)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dlauum");
    LAPACK_IMPL(dlauum)(UPLO,
                        N,
                        A,
                        LDA,
                        INFO);
}

//-- dopgtr --------------------------------------------------------------------
void
LAPACK_DECL(dopgtr)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    const DOUBLE     *TAU,
                    DOUBLE           *Q,
                    const INTEGER    *LDQ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dopgtr");
    LAPACK_IMPL(dopgtr)(UPLO,
                        N,
                        AP,
                        TAU,
                        Q,
                        LDQ,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dopmtr");
    LAPACK_IMPL(dopmtr)(SIDE,
                        UPLO,
                        TRANS,
                        M,
                        N,
                        AP,
                        TAU,
                        C,
                        LDC,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorbdb");
    LAPACK_IMPL(dorbdb)(TRANS,
                        SIGNS,
                        M,
                        P,
                        Q,
                        X11,
                        LDX11,
                        X12,
                        LDX12,
                        X21,
                        LDX21,
                        X22,
                        LDX22,
                        THETA,
                        PHI,
                        TAUP1,
                        TAUP2,
                        TAUQ1,
                        TAUQ2,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorcsd");
    LAPACK_IMPL(dorcsd)(JOBU1,
                        JOBU2,
                        JOBV1T,
                        JOBV2T,
                        TRANS,
                        SIGNS,
                        M,
                        P,
                        Q,
                        X11,
                        LDX11,
                        X12,
                        LDX12,
                        X21,
                        LDX21,
                        X22,
                        LDX22,
                        THETA,
                        U1,
                        LDU1,
                        U2,
                        LDU2,
                        V1T,
                        LDV1T,
                        V2T,
                        LDV2T,
                        WORK,
                        LWORK,
                        IWORK,
                        INFO);
}

//-- dorg2l --------------------------------------------------------------------
void
LAPACK_DECL(dorg2l)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorg2l");
    LAPACK_IMPL(dorg2l)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

//-- dorg2r --------------------------------------------------------------------
void
LAPACK_DECL(dorg2r)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorg2r");
    LAPACK_IMPL(dorg2r)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorgbr");
    LAPACK_IMPL(dorgbr)(VECT,
                        M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorghr");
    LAPACK_IMPL(dorghr)(N,
                        ILO,
                        IHI,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

//-- dorgl2 --------------------------------------------------------------------
void
LAPACK_DECL(dorgl2)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorgl2");
    LAPACK_IMPL(dorgl2)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

//-- dorglq --------------------------------------------------------------------
/*
void
LAPACK_DECL(dorglq)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorglq");
    LAPACK_IMPL(dorglq)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}
*/
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorgql");
    LAPACK_IMPL(dorgql)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

//-- dorgqr --------------------------------------------------------------------
/*
void
LAPACK_DECL(dorgqr)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorgqr");
    LAPACK_IMPL(dorgqr)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}
*/
//-- dorgr2 --------------------------------------------------------------------
void
LAPACK_DECL(dorgr2)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorgr2");
    LAPACK_IMPL(dorgr2)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorgrq");
    LAPACK_IMPL(dorgrq)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

//-- dorgtr --------------------------------------------------------------------
void
LAPACK_DECL(dorgtr)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorgtr");
    LAPACK_IMPL(dorgtr)(UPLO,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorm2l");
    LAPACK_IMPL(dorm2l)(SIDE,
                        TRANS,
                        M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        C,
                        LDC,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorm2r");
    LAPACK_IMPL(dorm2r)(SIDE,
                        TRANS,
                        M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        C,
                        LDC,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dormbr");
    LAPACK_IMPL(dormbr)(VECT,
                        SIDE,
                        TRANS,
                        M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        C,
                        LDC,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dormhr");
    LAPACK_IMPL(dormhr)(SIDE,
                        TRANS,
                        M,
                        N,
                        ILO,
                        IHI,
                        A,
                        LDA,
                        TAU,
                        C,
                        LDC,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dorml2");
    LAPACK_IMPL(dorml2)(SIDE,
                        TRANS,
                        M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        C,
                        LDC,
                        WORK,
                        INFO);
}

//-- dormlq --------------------------------------------------------------------
/*
void
LAPACK_DECL(dormlq)(const char       *SIDE,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dormlq");
    LAPACK_IMPL(dormlq)(SIDE,
                        TRANS,
                        M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        C,
                        LDC,
                        WORK,
                        LWORK,
                        INFO);
}
*/
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dormql");
    LAPACK_IMPL(dormql)(SIDE,
                        TRANS,
                        M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        C,
                        LDC,
                        WORK,
                        LWORK,
                        INFO);
}

//-- dormqr --------------------------------------------------------------------
/*
void
LAPACK_DECL(dormqr)(const char       *SIDE,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dormqr");
    LAPACK_IMPL(dormqr)(SIDE,
                        TRANS,
                        M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        C,
                        LDC,
                        WORK,
                        LWORK,
                        INFO);
}
*/
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dormr2");
    LAPACK_IMPL(dormr2)(SIDE,
                        TRANS,
                        M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        C,
                        LDC,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dormr3");
    LAPACK_IMPL(dormr3)(SIDE,
                        TRANS,
                        M,
                        N,
                        K,
                        L,
                        A,
                        LDA,
                        TAU,
                        C,
                        LDC,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dormrq");
    LAPACK_IMPL(dormrq)(SIDE,
                        TRANS,
                        M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        C,
                        LDC,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dormrz");
    LAPACK_IMPL(dormrz)(SIDE,
                        TRANS,
                        M,
                        N,
                        K,
                        L,
                        A,
                        LDA,
                        TAU,
                        C,
                        LDC,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dormtr");
    LAPACK_IMPL(dormtr)(SIDE,
                        UPLO,
                        TRANS,
                        M,
                        N,
                        A,
                        LDA,
                        TAU,
                        C,
                        LDC,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpbcon");
    LAPACK_IMPL(dpbcon)(UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        ANORM,
                        RCOND,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpbequ");
    LAPACK_IMPL(dpbequ)(UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        S,
                        SCOND,
                        AMAX,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpbrfs");
    LAPACK_IMPL(dpbrfs)(UPLO,
                        N,
                        KD,
                        NRHS,
                        AB,
                        LDAB,
                        AFB,
                        LDAFB,
                        B,
                        LDB,
                        X,
                        LDX,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dpbstf --------------------------------------------------------------------
void
LAPACK_DECL(dpbstf)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpbstf");
    LAPACK_IMPL(dpbstf)(UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dpbsv");
    LAPACK_IMPL(dpbsv)(UPLO,
                       N,
                       KD,
                       NRHS,
                       AB,
                       LDAB,
                       B,
                       LDB,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpbsvx");
    LAPACK_IMPL(dpbsvx)(FACT,
                        UPLO,
                        N,
                        KD,
                        NRHS,
                        AB,
                        LDAB,
                        AFB,
                        LDAFB,
                        EQUED,
                        S,
                        B,
                        LDB,
                        X,
                        LDX,
                        RCOND,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dpbtf2 --------------------------------------------------------------------
void
LAPACK_DECL(dpbtf2)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpbtf2");
    LAPACK_IMPL(dpbtf2)(UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        INFO);
}

//-- dpbtrf --------------------------------------------------------------------
void
LAPACK_DECL(dpbtrf)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE           *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpbtrf");
    LAPACK_IMPL(dpbtrf)(UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpbtrs");
    LAPACK_IMPL(dpbtrs)(UPLO,
                        N,
                        KD,
                        NRHS,
                        AB,
                        LDAB,
                        B,
                        LDB,
                        INFO);
}

//-- dpftrf --------------------------------------------------------------------
void
LAPACK_DECL(dpftrf)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpftrf");
    LAPACK_IMPL(dpftrf)(TRANSR,
                        UPLO,
                        N,
                        A,
                        INFO);
}

//-- dpftri --------------------------------------------------------------------
void
LAPACK_DECL(dpftri)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpftri");
    LAPACK_IMPL(dpftri)(TRANSR,
                        UPLO,
                        N,
                        A,
                        INFO);
}

//-- dpftrs --------------------------------------------------------------------
void
LAPACK_DECL(dpftrs)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpftrs");
    LAPACK_IMPL(dpftrs)(TRANSR,
                        UPLO,
                        N,
                        NRHS,
                        A,
                        B,
                        LDB,
                        INFO);
}

//-- dpocon --------------------------------------------------------------------
/*
void
LAPACK_DECL(dpocon)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpocon");
    LAPACK_IMPL(dpocon)(UPLO,
                        N,
                        A,
                        LDA,
                        ANORM,
                        RCOND,
                        WORK,
                        IWORK,
                        INFO);
}
*/
//-- dpoequ --------------------------------------------------------------------
void
LAPACK_DECL(dpoequ)(const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *S,
                    DOUBLE           *SCOND,
                    DOUBLE           *AMAX,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpoequ");
    LAPACK_IMPL(dpoequ)(N,
                        A,
                        LDA,
                        S,
                        SCOND,
                        AMAX,
                        INFO);
}

//-- dpoequb -------------------------------------------------------------------
void
LAPACK_DECL(dpoequb)(const INTEGER    *N,
                     const DOUBLE     *A,
                     const INTEGER    *LDA,
                     DOUBLE           *S,
                     DOUBLE           *SCOND,
                     DOUBLE           *AMAX,
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpoequb");
    LAPACK_IMPL(dpoequb)(N,
                         A,
                         LDA,
                         S,
                         SCOND,
                         AMAX,
                         INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dporfs");
    LAPACK_IMPL(dporfs)(UPLO,
                        N,
                        NRHS,
                        A,
                        LDA,
                        AF,
                        LDAF,
                        B,
                        LDB,
                        X,
                        LDX,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dposv ---------------------------------------------------------------------
/*
void
LAPACK_DECL(dposv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dposv");
    LAPACK_IMPL(dposv)(UPLO,
                       N,
                       NRHS,
                       A,
                       LDA,
                       B,
                       LDB,
                       INFO);
}
*/
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dposvx");
    LAPACK_IMPL(dposvx)(FACT,
                        UPLO,
                        N,
                        NRHS,
                        A,
                        LDA,
                        AF,
                        LDAF,
                        EQUED,
                        S,
                        B,
                        LDB,
                        X,
                        LDX,
                        RCOND,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dpotf2 --------------------------------------------------------------------
void
LAPACK_DECL(dpotf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpotf2");
    LAPACK_IMPL(dpotf2)(UPLO,
                        N,
                        A,
                        LDA,
                        INFO);
}

//-- dpotrf --------------------------------------------------------------------
/*
void
LAPACK_DECL(dpotrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpotrf");
    LAPACK_IMPL(dpotrf)(UPLO,
                        N,
                        A,
                        LDA,
                        INFO);
}
*/
//-- dpotri --------------------------------------------------------------------
/*
void
LAPACK_DECL(dpotri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpotri");
    LAPACK_IMPL(dpotri)(UPLO,
                        N,
                        A,
                        LDA,
                        INFO);
}
*/
//-- dpotrs --------------------------------------------------------------------
/*
void
LAPACK_DECL(dpotrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpotrs");
    LAPACK_IMPL(dpotrs)(UPLO,
                        N,
                        NRHS,
                        A,
                        LDA,
                        B,
                        LDB,
                        INFO);
}
*/
//-- dppcon --------------------------------------------------------------------
void
LAPACK_DECL(dppcon)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dppcon");
    LAPACK_IMPL(dppcon)(UPLO,
                        N,
                        AP,
                        ANORM,
                        RCOND,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dppequ --------------------------------------------------------------------
void
LAPACK_DECL(dppequ)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *S,
                    DOUBLE           *SCOND,
                    DOUBLE           *AMAX,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dppequ");
    LAPACK_IMPL(dppequ)(UPLO,
                        N,
                        AP,
                        S,
                        SCOND,
                        AMAX,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpprfs");
    LAPACK_IMPL(dpprfs)(UPLO,
                        N,
                        NRHS,
                        AP,
                        AFP,
                        B,
                        LDB,
                        X,
                        LDX,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dppsv ---------------------------------------------------------------------
void
LAPACK_DECL(dppsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *AP,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dppsv");
    LAPACK_IMPL(dppsv)(UPLO,
                       N,
                       NRHS,
                       AP,
                       B,
                       LDB,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dppsvx");
    LAPACK_IMPL(dppsvx)(FACT,
                        UPLO,
                        N,
                        NRHS,
                        AP,
                        AFP,
                        EQUED,
                        S,
                        B,
                        LDB,
                        X,
                        LDX,
                        RCOND,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dpptrf --------------------------------------------------------------------
void
LAPACK_DECL(dpptrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpptrf");
    LAPACK_IMPL(dpptrf)(UPLO,
                        N,
                        AP,
                        INFO);
}

//-- dpptri --------------------------------------------------------------------
void
LAPACK_DECL(dpptri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpptri");
    LAPACK_IMPL(dpptri)(UPLO,
                        N,
                        AP,
                        INFO);
}

//-- dpptrs --------------------------------------------------------------------
void
LAPACK_DECL(dpptrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AP,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpptrs");
    LAPACK_IMPL(dpptrs)(UPLO,
                        N,
                        NRHS,
                        AP,
                        B,
                        LDB,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpstf2");
    LAPACK_IMPL(dpstf2)(UPLO,
                        N,
                        A,
                        LDA,
                        PIV,
                        RANK,
                        TOL,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpstrf");
    LAPACK_IMPL(dpstrf)(UPLO,
                        N,
                        A,
                        LDA,
                        PIV,
                        RANK,
                        TOL,
                        WORK,
                        INFO);
}

//-- dptcon --------------------------------------------------------------------
void
LAPACK_DECL(dptcon)(const INTEGER    *N,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dptcon");
    LAPACK_IMPL(dptcon)(N,
                        D,
                        E,
                        ANORM,
                        RCOND,
                        WORK,
                        INFO);
}

//-- dpteqr --------------------------------------------------------------------
void
LAPACK_DECL(dpteqr)(const char       *COMPZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpteqr");
    LAPACK_IMPL(dpteqr)(COMPZ,
                        N,
                        D,
                        E,
                        Z,
                        LDZ,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dptrfs");
    LAPACK_IMPL(dptrfs)(N,
                        NRHS,
                        D,
                        E,
                        DF,
                        EF,
                        B,
                        LDB,
                        X,
                        LDX,
                        FERR,
                        BERR,
                        WORK,
                        INFO);
}

//-- dptsv ---------------------------------------------------------------------
void
LAPACK_DECL(dptsv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *D,
                   DOUBLE               *E,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dptsv");
    LAPACK_IMPL(dptsv)(N,
                       NRHS,
                       D,
                       E,
                       B,
                       LDB,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dptsvx");
    LAPACK_IMPL(dptsvx)(FACT,
                        N,
                        NRHS,
                        D,
                        E,
                        DF,
                        EF,
                        B,
                        LDB,
                        X,
                        LDX,
                        RCOND,
                        FERR,
                        BERR,
                        WORK,
                        INFO);
}

//-- dpttrf --------------------------------------------------------------------
void
LAPACK_DECL(dpttrf)(const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpttrf");
    LAPACK_IMPL(dpttrf)(N,
                        D,
                        E,
                        INFO);
}

//-- dpttrs --------------------------------------------------------------------
void
LAPACK_DECL(dpttrs)(const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dpttrs");
    LAPACK_IMPL(dpttrs)(N,
                        NRHS,
                        D,
                        E,
                        B,
                        LDB,
                        INFO);
}

//-- dptts2 --------------------------------------------------------------------
void
LAPACK_DECL(dptts2)(const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *D,
                    const DOUBLE     *E,
                    DOUBLE           *B,
                    const INTEGER    *LDB)
{
    DEBUG_LAPACK_STUB("dptts2");
    LAPACK_IMPL(dptts2)(N,
                        NRHS,
                        D,
                        E,
                        B,
                        LDB);
}

//-- drscl ---------------------------------------------------------------------
void
LAPACK_DECL(drscl)(const INTEGER        *N,
                   const DOUBLE         *SA,
                   DOUBLE               *SX,
                   const INTEGER        *INCX)
{
    DEBUG_LAPACK_STUB("drscl");
    LAPACK_IMPL(drscl)(N,
                       SA,
                       SX,
                       INCX);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dsbev");
    LAPACK_IMPL(dsbev)(JOBZ,
                       UPLO,
                       N,
                       KD,
                       AB,
                       LDAB,
                       W,
                       Z,
                       LDZ,
                       WORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsbevd");
    LAPACK_IMPL(dsbevd)(JOBZ,
                        UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        W,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsbevx");
    LAPACK_IMPL(dsbevx)(JOBZ,
                        RANGE,
                        UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        Q,
                        LDQ,
                        VL,
                        VU,
                        IL,
                        IU,
                        ABSTOL,
                        M,
                        W,
                        Z,
                        LDZ,
                        WORK,
                        IWORK,
                        IFAIL,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsbgst");
    LAPACK_IMPL(dsbgst)(VECT,
                        UPLO,
                        N,
                        KA,
                        KB,
                        AB,
                        LDAB,
                        BB,
                        LDBB,
                        X,
                        LDX,
                        WORK,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dsbgv");
    LAPACK_IMPL(dsbgv)(JOBZ,
                       UPLO,
                       N,
                       KA,
                       KB,
                       AB,
                       LDAB,
                       BB,
                       LDBB,
                       W,
                       Z,
                       LDZ,
                       WORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsbgvd");
    LAPACK_IMPL(dsbgvd)(JOBZ,
                        UPLO,
                        N,
                        KA,
                        KB,
                        AB,
                        LDAB,
                        BB,
                        LDBB,
                        W,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsbgvx");
    LAPACK_IMPL(dsbgvx)(JOBZ,
                        RANGE,
                        UPLO,
                        N,
                        KA,
                        KB,
                        AB,
                        LDAB,
                        BB,
                        LDBB,
                        Q,
                        LDQ,
                        VL,
                        VU,
                        IL,
                        IU,
                        ABSTOL,
                        M,
                        W,
                        Z,
                        LDZ,
                        WORK,
                        IWORK,
                        IFAIL,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsbtrd");
    LAPACK_IMPL(dsbtrd)(VECT,
                        UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        D,
                        E,
                        Q,
                        LDQ,
                        WORK,
                        INFO);
}

//-- dsecnd --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dsecnd)()
{
    DEBUG_LAPACK_STUB("dsecnd");
    return LAPACK_IMPL(dsecnd)();
}

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
                   DOUBLE               *C)
{
    DEBUG_LAPACK_STUB("dsfrk");
    LAPACK_IMPL(dsfrk)(TRANSR,
                       UPLO,
                       TRANS,
                       N,
                       K,
                       ALPHA,
                       A,
                       LDA,
                       BETA,
                       C);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsgesv");
    LAPACK_IMPL(dsgesv)(N,
                        NRHS,
                        A,
                        LDA,
                        IPIV,
                        B,
                        LDB,
                        X,
                        LDX,
                        WORK,
                        SWORK,
                        ITER,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dspcon");
    LAPACK_IMPL(dspcon)(UPLO,
                        N,
                        AP,
                        IPIV,
                        ANORM,
                        RCOND,
                        WORK,
                        IWORK,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dspev");
    LAPACK_IMPL(dspev)(JOBZ,
                       UPLO,
                       N,
                       AP,
                       W,
                       Z,
                       LDZ,
                       WORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dspevd");
    LAPACK_IMPL(dspevd)(JOBZ,
                        UPLO,
                        N,
                        AP,
                        W,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dspevx");
    LAPACK_IMPL(dspevx)(JOBZ,
                        RANGE,
                        UPLO,
                        N,
                        AP,
                        VL,
                        VU,
                        IL,
                        IU,
                        ABSTOL,
                        M,
                        W,
                        Z,
                        LDZ,
                        WORK,
                        IWORK,
                        IFAIL,
                        INFO);
}

//-- dspgst --------------------------------------------------------------------
void
LAPACK_DECL(dspgst)(const INTEGER    *ITYPE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    const DOUBLE     *BP,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dspgst");
    LAPACK_IMPL(dspgst)(ITYPE,
                        UPLO,
                        N,
                        AP,
                        BP,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dspgv");
    LAPACK_IMPL(dspgv)(ITYPE,
                       JOBZ,
                       UPLO,
                       N,
                       AP,
                       BP,
                       W,
                       Z,
                       LDZ,
                       WORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dspgvd");
    LAPACK_IMPL(dspgvd)(ITYPE,
                        JOBZ,
                        UPLO,
                        N,
                        AP,
                        BP,
                        W,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dspgvx");
    LAPACK_IMPL(dspgvx)(ITYPE,
                        JOBZ,
                        RANGE,
                        UPLO,
                        N,
                        AP,
                        BP,
                        VL,
                        VU,
                        IL,
                        IU,
                        ABSTOL,
                        M,
                        W,
                        Z,
                        LDZ,
                        WORK,
                        IWORK,
                        IFAIL,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsposv");
    LAPACK_IMPL(dsposv)(UPLO,
                        N,
                        NRHS,
                        A,
                        LDA,
                        B,
                        LDB,
                        X,
                        LDX,
                        WORK,
                        SWORK,
                        ITER,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsprfs");
    LAPACK_IMPL(dsprfs)(UPLO,
                        N,
                        NRHS,
                        AP,
                        AFP,
                        IPIV,
                        B,
                        LDB,
                        X,
                        LDX,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dspsv ---------------------------------------------------------------------
void
LAPACK_DECL(dspsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *AP,
                   INTEGER              *IPIV,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dspsv");
    LAPACK_IMPL(dspsv)(UPLO,
                       N,
                       NRHS,
                       AP,
                       IPIV,
                       B,
                       LDB,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dspsvx");
    LAPACK_IMPL(dspsvx)(FACT,
                        UPLO,
                        N,
                        NRHS,
                        AP,
                        AFP,
                        IPIV,
                        B,
                        LDB,
                        X,
                        LDX,
                        RCOND,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dsptrd --------------------------------------------------------------------
void
LAPACK_DECL(dsptrd)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *TAU,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsptrd");
    LAPACK_IMPL(dsptrd)(UPLO,
                        N,
                        AP,
                        D,
                        E,
                        TAU,
                        INFO);
}

//-- dsptrf --------------------------------------------------------------------
void
LAPACK_DECL(dsptrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsptrf");
    LAPACK_IMPL(dsptrf)(UPLO,
                        N,
                        AP,
                        IPIV,
                        INFO);
}

//-- dsptri --------------------------------------------------------------------
void
LAPACK_DECL(dsptri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    const INTEGER    *IPIV,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsptri");
    LAPACK_IMPL(dsptri)(UPLO,
                        N,
                        AP,
                        IPIV,
                        WORK,
                        INFO);
}

//-- dsptrs --------------------------------------------------------------------
void
LAPACK_DECL(dsptrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *AP,
                    const INTEGER    *IPIV,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsptrs");
    LAPACK_IMPL(dsptrs)(UPLO,
                        N,
                        NRHS,
                        AP,
                        IPIV,
                        B,
                        LDB,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dstebz");
    LAPACK_IMPL(dstebz)(RANGE,
                        ORDER,
                        N,
                        VL,
                        VU,
                        IL,
                        IU,
                        ABSTOL,
                        D,
                        E,
                        M,
                        NSPLIT,
                        W,
                        IBLOCK,
                        ISPLIT,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dstedc");
    LAPACK_IMPL(dstedc)(COMPZ,
                        N,
                        D,
                        E,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dstegr");
    LAPACK_IMPL(dstegr)(JOBZ,
                        RANGE,
                        N,
                        D,
                        E,
                        VL,
                        VU,
                        IL,
                        IU,
                        ABSTOL,
                        M,
                        W,
                        Z,
                        LDZ,
                        ISUPPZ,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dstein");
    LAPACK_IMPL(dstein)(N,
                        D,
                        E,
                        M,
                        W,
                        IBLOCK,
                        ISPLIT,
                        Z,
                        LDZ,
                        WORK,
                        IWORK,
                        IFAIL,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dstemr");
    LAPACK_IMPL(dstemr)(JOBZ,
                        RANGE,
                        N,
                        D,
                        E,
                        VL,
                        VU,
                        IL,
                        IU,
                        M,
                        W,
                        Z,
                        LDZ,
                        NZC,
                        ISUPPZ,
                        TRYRAC,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

//-- dsteqr --------------------------------------------------------------------
void
LAPACK_DECL(dsteqr)(const char       *COMPZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsteqr");
    LAPACK_IMPL(dsteqr)(COMPZ,
                        N,
                        D,
                        E,
                        Z,
                        LDZ,
                        WORK,
                        INFO);
}

//-- dsterf --------------------------------------------------------------------
void
LAPACK_DECL(dsterf)(const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsterf");
    LAPACK_IMPL(dsterf)(N,
                        D,
                        E,
                        INFO);
}

//-- dstev ---------------------------------------------------------------------
void
LAPACK_DECL(dstev)(const char           *JOBZ,
                   const INTEGER        *N,
                   DOUBLE               *D,
                   DOUBLE               *E,
                   DOUBLE               *Z,
                   const INTEGER        *LDZ,
                   DOUBLE               *WORK,
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dstev");
    LAPACK_IMPL(dstev)(JOBZ,
                       N,
                       D,
                       E,
                       Z,
                       LDZ,
                       WORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dstevd");
    LAPACK_IMPL(dstevd)(JOBZ,
                        N,
                        D,
                        E,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dstevr");
    LAPACK_IMPL(dstevr)(JOBZ,
                        RANGE,
                        N,
                        D,
                        E,
                        VL,
                        VU,
                        IL,
                        IU,
                        ABSTOL,
                        M,
                        W,
                        Z,
                        LDZ,
                        ISUPPZ,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dstevx");
    LAPACK_IMPL(dstevx)(JOBZ,
                        RANGE,
                        N,
                        D,
                        E,
                        VL,
                        VU,
                        IL,
                        IU,
                        ABSTOL,
                        M,
                        W,
                        Z,
                        LDZ,
                        WORK,
                        IWORK,
                        IFAIL,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsycon");
    LAPACK_IMPL(dsycon)(UPLO,
                        N,
                        A,
                        LDA,
                        IPIV,
                        ANORM,
                        RCOND,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dsyconv -------------------------------------------------------------------
void
LAPACK_DECL(dsyconv)(const char       *UPLO,
                     const char       *WAY,
                     const INTEGER    *N,
                     const DOUBLE     *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE           *WORK,
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsyconv");
    LAPACK_IMPL(dsyconv)(UPLO,
                         WAY,
                         N,
                         A,
                         LDA,
                         IPIV,
                         WORK,
                         INFO);
}

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
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsyequb");
    LAPACK_IMPL(dsyequb)(UPLO,
                         N,
                         A,
                         LDA,
                         S,
                         SCOND,
                         AMAX,
                         WORK,
                         INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dsyev");
    LAPACK_IMPL(dsyev)(JOBZ,
                       UPLO,
                       N,
                       A,
                       LDA,
                       W,
                       WORK,
                       LWORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsyevd");
    LAPACK_IMPL(dsyevd)(JOBZ,
                        UPLO,
                        N,
                        A,
                        LDA,
                        W,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsyevr");
    LAPACK_IMPL(dsyevr)(JOBZ,
                        RANGE,
                        UPLO,
                        N,
                        A,
                        LDA,
                        VL,
                        VU,
                        IL,
                        IU,
                        ABSTOL,
                        M,
                        W,
                        Z,
                        LDZ,
                        ISUPPZ,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsyevx");
    LAPACK_IMPL(dsyevx)(JOBZ,
                        RANGE,
                        UPLO,
                        N,
                        A,
                        LDA,
                        VL,
                        VU,
                        IL,
                        IU,
                        ABSTOL,
                        M,
                        W,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        IWORK,
                        IFAIL,
                        INFO);
}

//-- dsygs2 --------------------------------------------------------------------
void
LAPACK_DECL(dsygs2)(const INTEGER    *ITYPE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsygs2");
    LAPACK_IMPL(dsygs2)(ITYPE,
                        UPLO,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        INFO);
}

//-- dsygst --------------------------------------------------------------------
void
LAPACK_DECL(dsygst)(const INTEGER    *ITYPE,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsygst");
    LAPACK_IMPL(dsygst)(ITYPE,
                        UPLO,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dsygv");
    LAPACK_IMPL(dsygv)(ITYPE,
                       JOBZ,
                       UPLO,
                       N,
                       A,
                       LDA,
                       B,
                       LDB,
                       W,
                       WORK,
                       LWORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsygvd");
    LAPACK_IMPL(dsygvd)(ITYPE,
                        JOBZ,
                        UPLO,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        W,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsygvx");
    LAPACK_IMPL(dsygvx)(ITYPE,
                        JOBZ,
                        RANGE,
                        UPLO,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        VL,
                        VU,
                        IL,
                        IU,
                        ABSTOL,
                        M,
                        W,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        IWORK,
                        IFAIL,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsyrfs");
    LAPACK_IMPL(dsyrfs)(UPLO,
                        N,
                        NRHS,
                        A,
                        LDA,
                        AF,
                        LDAF,
                        IPIV,
                        B,
                        LDB,
                        X,
                        LDX,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("dsysv");
    LAPACK_IMPL(dsysv)(UPLO,
                       N,
                       NRHS,
                       A,
                       LDA,
                       IPIV,
                       B,
                       LDB,
                       WORK,
                       LWORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsysvx");
    LAPACK_IMPL(dsysvx)(FACT,
                        UPLO,
                        N,
                        NRHS,
                        A,
                        LDA,
                        AF,
                        LDAF,
                        IPIV,
                        B,
                        LDB,
                        X,
                        LDX,
                        RCOND,
                        FERR,
                        BERR,
                        WORK,
                        LWORK,
                        IWORK,
                        INFO);
}

//-- dsyswapr ------------------------------------------------------------------
void
LAPACK_DECL(dsyswapr)(const char       *UPLO,
                      const INTEGER    *N,
                      DOUBLE           *A,
                      const INTEGER    *LDA,
                      const INTEGER    *I1,
                      const INTEGER    *I2)
{
    DEBUG_LAPACK_STUB("dsyswapr");
    LAPACK_IMPL(dsyswapr)(UPLO,
                          N,
                          A,
                          LDA,
                          I1,
                          I2);
}

//-- dsytd2 --------------------------------------------------------------------
void
LAPACK_DECL(dsytd2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE           *TAU,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsytd2");
    LAPACK_IMPL(dsytd2)(UPLO,
                        N,
                        A,
                        LDA,
                        D,
                        E,
                        TAU,
                        INFO);
}

//-- dsytf2 --------------------------------------------------------------------
void
LAPACK_DECL(dsytf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsytf2");
    LAPACK_IMPL(dsytf2)(UPLO,
                        N,
                        A,
                        LDA,
                        IPIV,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsytrd");
    LAPACK_IMPL(dsytrd)(UPLO,
                        N,
                        A,
                        LDA,
                        D,
                        E,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

//-- dsytrf --------------------------------------------------------------------
void
LAPACK_DECL(dsytrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsytrf");
    LAPACK_IMPL(dsytrf)(UPLO,
                        N,
                        A,
                        LDA,
                        IPIV,
                        WORK,
                        LWORK,
                        INFO);
}

//-- dsytri --------------------------------------------------------------------
void
LAPACK_DECL(dsytri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsytri");
    LAPACK_IMPL(dsytri)(UPLO,
                        N,
                        A,
                        LDA,
                        IPIV,
                        WORK,
                        INFO);
}

//-- dsytri2 -------------------------------------------------------------------
void
LAPACK_DECL(dsytri2)(const char       *UPLO,
                     const INTEGER    *N,
                     DOUBLE           *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE           *WORK,
                     const INTEGER    *LWORK,
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsytri2");
    LAPACK_IMPL(dsytri2)(UPLO,
                         N,
                         A,
                         LDA,
                         IPIV,
                         WORK,
                         LWORK,
                         INFO);
}

//-- dsytri2x ------------------------------------------------------------------
void
LAPACK_DECL(dsytri2x)(const char       *UPLO,
                      const INTEGER    *N,
                      DOUBLE           *A,
                      const INTEGER    *LDA,
                      const INTEGER    *IPIV,
                      DOUBLE           *WORK,
                      const INTEGER    *NB,
                      INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsytri2x");
    LAPACK_IMPL(dsytri2x)(UPLO,
                          N,
                          A,
                          LDA,
                          IPIV,
                          WORK,
                          NB,
                          INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsytrs");
    LAPACK_IMPL(dsytrs)(UPLO,
                        N,
                        NRHS,
                        A,
                        LDA,
                        IPIV,
                        B,
                        LDB,
                        INFO);
}

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
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dsytrs2");
    LAPACK_IMPL(dsytrs2)(UPLO,
                         N,
                         NRHS,
                         A,
                         LDA,
                         IPIV,
                         B,
                         LDB,
                         WORK,
                         INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtbcon");
    LAPACK_IMPL(dtbcon)(NORM,
                        UPLO,
                        DIAG,
                        N,
                        KD,
                        AB,
                        LDAB,
                        RCOND,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtbrfs");
    LAPACK_IMPL(dtbrfs)(UPLO,
                        TRANS,
                        DIAG,
                        N,
                        KD,
                        NRHS,
                        AB,
                        LDAB,
                        B,
                        LDB,
                        X,
                        LDX,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtbtrs");
    LAPACK_IMPL(dtbtrs)(UPLO,
                        TRANS,
                        DIAG,
                        N,
                        KD,
                        NRHS,
                        AB,
                        LDAB,
                        B,
                        LDB,
                        INFO);
}

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
                   const INTEGER        *LDB)
{
    DEBUG_LAPACK_STUB("dtfsm");
    LAPACK_IMPL(dtfsm)(TRANSR,
                       SIDE,
                       UPLO,
                       TRANS,
                       DIAG,
                       M,
                       N,
                       ALPHA,
                       A,
                       B,
                       LDB);
}

//-- dtftri --------------------------------------------------------------------
void
LAPACK_DECL(dtftri)(const char       *TRANSR,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtftri");
    LAPACK_IMPL(dtftri)(TRANSR,
                        UPLO,
                        DIAG,
                        N,
                        A,
                        INFO);
}

//-- dtfttp --------------------------------------------------------------------
void
LAPACK_DECL(dtfttp)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *ARF,
                    DOUBLE           *AP,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtfttp");
    LAPACK_IMPL(dtfttp)(TRANSR,
                        UPLO,
                        N,
                        ARF,
                        AP,
                        INFO);
}

//-- dtfttr --------------------------------------------------------------------
void
LAPACK_DECL(dtfttr)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *ARF,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtfttr");
    LAPACK_IMPL(dtfttr)(TRANSR,
                        UPLO,
                        N,
                        ARF,
                        A,
                        LDA,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtgevc");
    LAPACK_IMPL(dtgevc)(SIDE,
                        HOWMNY,
                        SELECT,
                        N,
                        S,
                        LDS,
                        P,
                        LDP,
                        VL,
                        LDVL,
                        VR,
                        LDVR,
                        MM,
                        M,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtgex2");
    LAPACK_IMPL(dtgex2)(WANTQ,
                        WANTZ,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        Q,
                        LDQ,
                        Z,
                        LDZ,
                        J1,
                        N1,
                        N2,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtgexc");
    LAPACK_IMPL(dtgexc)(WANTQ,
                        WANTZ,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        Q,
                        LDQ,
                        Z,
                        LDZ,
                        IFST,
                        ILST,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtgsen");
    LAPACK_IMPL(dtgsen)(IJOB,
                        WANTQ,
                        WANTZ,
                        SELECT,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        ALPHAR,
                        ALPHAI,
                        BETA,
                        Q,
                        LDQ,
                        Z,
                        LDZ,
                        M,
                        PL,
                        PR,
                        DIF,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtgsja");
    LAPACK_IMPL(dtgsja)(JOBU,
                        JOBV,
                        JOBQ,
                        M,
                        P,
                        N,
                        K,
                        L,
                        A,
                        LDA,
                        B,
                        LDB,
                        TOLA,
                        TOLB,
                        ALPHA,
                        BETA,
                        U,
                        LDU,
                        V,
                        LDV,
                        Q,
                        LDQ,
                        WORK,
                        NCYCLE,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtgsna");
    LAPACK_IMPL(dtgsna)(JOB,
                        HOWMNY,
                        SELECT,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        VL,
                        LDVL,
                        VR,
                        LDVR,
                        S,
                        DIF,
                        MM,
                        M,
                        WORK,
                        LWORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtgsy2");
    LAPACK_IMPL(dtgsy2)(TRANS,
                        IJOB,
                        M,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        C,
                        LDC,
                        D,
                        LDD,
                        E,
                        LDE,
                        F,
                        LDF,
                        SCALE,
                        RDSUM,
                        RDSCAL,
                        IWORK,
                        PQ,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtgsyl");
    LAPACK_IMPL(dtgsyl)(TRANS,
                        IJOB,
                        M,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        C,
                        LDC,
                        D,
                        LDD,
                        E,
                        LDE,
                        F,
                        LDF,
                        SCALE,
                        DIF,
                        WORK,
                        LWORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtpcon");
    LAPACK_IMPL(dtpcon)(NORM,
                        UPLO,
                        DIAG,
                        N,
                        AP,
                        RCOND,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtprfs");
    LAPACK_IMPL(dtprfs)(UPLO,
                        TRANS,
                        DIAG,
                        N,
                        NRHS,
                        AP,
                        B,
                        LDB,
                        X,
                        LDX,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

//-- dtptri --------------------------------------------------------------------
void
LAPACK_DECL(dtptri)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE           *AP,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtptri");
    LAPACK_IMPL(dtptri)(UPLO,
                        DIAG,
                        N,
                        AP,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtptrs");
    LAPACK_IMPL(dtptrs)(UPLO,
                        TRANS,
                        DIAG,
                        N,
                        NRHS,
                        AP,
                        B,
                        LDB,
                        INFO);
}

//-- dtpttf --------------------------------------------------------------------
void
LAPACK_DECL(dtpttf)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *ARF,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtpttf");
    LAPACK_IMPL(dtpttf)(TRANSR,
                        UPLO,
                        N,
                        AP,
                        ARF,
                        INFO);
}

//-- dtpttr --------------------------------------------------------------------
void
LAPACK_DECL(dtpttr)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *AP,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtpttr");
    LAPACK_IMPL(dtpttr)(UPLO,
                        N,
                        AP,
                        A,
                        LDA,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtrcon");
    LAPACK_IMPL(dtrcon)(NORM,
                        UPLO,
                        DIAG,
                        N,
                        A,
                        LDA,
                        RCOND,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtrevc");
    LAPACK_IMPL(dtrevc)(SIDE,
                        HOWMNY,
                        SELECT,
                        N,
                        T,
                        LDT,
                        VL,
                        LDVL,
                        VR,
                        LDVR,
                        MM,
                        M,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtrexc");
    LAPACK_IMPL(dtrexc)(COMPQ,
                        N,
                        T,
                        LDT,
                        Q,
                        LDQ,
                        IFST,
                        ILST,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtrrfs");
    LAPACK_IMPL(dtrrfs)(UPLO,
                        TRANS,
                        DIAG,
                        N,
                        NRHS,
                        A,
                        LDA,
                        B,
                        LDB,
                        X,
                        LDX,
                        FERR,
                        BERR,
                        WORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtrsen");
    LAPACK_IMPL(dtrsen)(JOB,
                        COMPQ,
                        SELECT,
                        N,
                        T,
                        LDT,
                        Q,
                        LDQ,
                        WR,
                        WI,
                        M,
                        S,
                        SEP,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtrsna");
    LAPACK_IMPL(dtrsna)(JOB,
                        HOWMNY,
                        SELECT,
                        N,
                        T,
                        LDT,
                        VL,
                        LDVL,
                        VR,
                        LDVR,
                        S,
                        SEP,
                        MM,
                        M,
                        WORK,
                        LDWORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtrsyl");
    LAPACK_IMPL(dtrsyl)(TRANA,
                        TRANB,
                        ISGN,
                        M,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        C,
                        LDC,
                        SCALE,
                        INFO);
}

//-- dtrti2 --------------------------------------------------------------------
void
LAPACK_DECL(dtrti2)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtrti2");
    LAPACK_IMPL(dtrti2)(UPLO,
                        DIAG,
                        N,
                        A,
                        LDA,
                        INFO);
}

//-- dtrtri --------------------------------------------------------------------
/*
void
LAPACK_DECL(dtrtri)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtrtri");
    LAPACK_IMPL(dtrtri)(UPLO,
                        DIAG,
                        N,
                        A,
                        LDA,
                        INFO);
}
*/
//-- dtrtrs --------------------------------------------------------------------
/*
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtrtrs");
    LAPACK_IMPL(dtrtrs)(UPLO,
                        TRANS,
                        DIAG,
                        N,
                        NRHS,
                        A,
                        LDA,
                        B,
                        LDB,
                        INFO);
}
*/
//-- dtrttf --------------------------------------------------------------------
void
LAPACK_DECL(dtrttf)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *ARF,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtrttf");
    LAPACK_IMPL(dtrttf)(TRANSR,
                        UPLO,
                        N,
                        A,
                        LDA,
                        ARF,
                        INFO);
}

//-- dtrttp --------------------------------------------------------------------
void
LAPACK_DECL(dtrttp)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *AP,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtrttp");
    LAPACK_IMPL(dtrttp)(UPLO,
                        N,
                        A,
                        LDA,
                        AP,
                        INFO);
}

//-- dtzrqf --------------------------------------------------------------------
void
LAPACK_DECL(dtzrqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtzrqf");
    LAPACK_IMPL(dtzrqf)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        INFO);
}

//-- dtzrzf --------------------------------------------------------------------
void
LAPACK_DECL(dtzrzf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("dtzrzf");
    LAPACK_IMPL(dtzrzf)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

//-- dzsum1 --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(dzsum1)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *CX,
                    const INTEGER            *INCX)
{
    DEBUG_LAPACK_STUB("dzsum1");
    return LAPACK_IMPL(dzsum1)(N,
                               CX,
                               INCX);
}

//-- ieeeck --------------------------------------------------------------------
INTEGER
LAPACK_DECL(ieeeck)(const INTEGER    *ISPEC,
                    const FLOAT      *ZERO,
                    const FLOAT      *ONE)
{
    DEBUG_LAPACK_STUB("ieeeck");
    return LAPACK_IMPL(ieeeck)(ISPEC,
                               ZERO,
                               ONE);
}

//-- iladlc --------------------------------------------------------------------
INTEGER
LAPACK_DECL(iladlc)(const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA)
{
    DEBUG_LAPACK_STUB("iladlc");
    return LAPACK_IMPL(iladlc)(M,
                               N,
                               A,
                               LDA);
}

//-- iladlr --------------------------------------------------------------------
INTEGER
LAPACK_DECL(iladlr)(const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA)
{
    DEBUG_LAPACK_STUB("iladlr");
    return LAPACK_IMPL(iladlr)(M,
                               N,
                               A,
                               LDA);
}

//-- ilaslc --------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilaslc)(const INTEGER    *M,
                    const INTEGER    *N,
                    const FLOAT      *A,
                    const INTEGER    *LDA)
{
    DEBUG_LAPACK_STUB("ilaslc");
    return LAPACK_IMPL(ilaslc)(M,
                               N,
                               A,
                               LDA);
}

//-- ilaslr --------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilaslr)(const INTEGER    *M,
                    const INTEGER    *N,
                    const FLOAT      *A,
                    const INTEGER    *LDA)
{
    DEBUG_LAPACK_STUB("ilaslr");
    return LAPACK_IMPL(ilaslr)(M,
                               N,
                               A,
                               LDA);
}

//-- ilatrans ------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilatrans)(const char       *TRANS)
{
    DEBUG_LAPACK_STUB("ilatrans");
    return LAPACK_IMPL(ilatrans)(TRANS);
}

//-- ilauplo -------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilauplo)(const char   *UPLO)
{
    DEBUG_LAPACK_STUB("ilauplo");
    return LAPACK_IMPL(ilauplo)(UPLO);
}

//-- ilaver --------------------------------------------------------------------
void
LAPACK_DECL(ilaver)(INTEGER  *VERS_MAJOR,
                    INTEGER  *VERS_MINOR,
                    INTEGER  *VERS_PATCH)
{
    DEBUG_LAPACK_STUB("ilaver");
    LAPACK_IMPL(ilaver)(VERS_MAJOR,
                        VERS_MINOR,
                        VERS_PATCH);
}

//-- lsame ---------------------------------------------------------------------
/*
LOGICAL
LAPACK_DECL(lsame)(const char       *CA,
                   const char       *CB)
{
    DEBUG_LAPACK_STUB("lsame");
    return LAPACK_IMPL(lsame)(CA,
                              CB);
}
*/
//-- lsamen --------------------------------------------------------------------
/*
LOGICAL
LAPACK_DECL(lsamen)(const INTEGER    *N,
                    const char       *CA,
                    const char       *CB)
{
    DEBUG_LAPACK_STUB("lsamen");
    return LAPACK_IMPL(lsamen)(N,
                               CA,
                               CB);
}
*/
//-- sgetf2 --------------------------------------------------------------------
void
LAPACK_DECL(sgetf2)(const INTEGER    *M,
                    const INTEGER    *N,
                    FLOAT            *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("sgetf2");
    LAPACK_IMPL(sgetf2)(M,
                        N,
                        A,
                        LDA,
                        IPIV,
                        INFO);
}

//-- sgetrf --------------------------------------------------------------------
void
LAPACK_DECL(sgetrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    FLOAT            *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("sgetrf");
    LAPACK_IMPL(sgetrf)(M,
                        N,
                        A,
                        LDA,
                        IPIV,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("sgetrs");
    LAPACK_IMPL(sgetrs)(TRANS,
                        N,
                        NRHS,
                        A,
                        LDA,
                        IPIV,
                        B,
                        LDB,
                        INFO);
}

//-- sisnan --------------------------------------------------------------------
LOGICAL
LAPACK_DECL(sisnan)(const FLOAT  *SIN)
{
    DEBUG_LAPACK_STUB("sisnan");
    return LAPACK_IMPL(sisnan)(SIN);
}

//-- slag2d --------------------------------------------------------------------
void
LAPACK_DECL(slag2d)(const INTEGER    *M,
                    const INTEGER    *N,
                    const FLOAT      *SA,
                    const INTEGER    *LDSA,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("slag2d");
    LAPACK_IMPL(slag2d)(M,
                        N,
                        SA,
                        LDSA,
                        A,
                        LDA,
                        INFO);
}

//-- slaisnan ------------------------------------------------------------------
LOGICAL
LAPACK_DECL(slaisnan)(const FLOAT      *SIN1,
                      const FLOAT      *SIN2)
{
    DEBUG_LAPACK_STUB("slaisnan");
    return LAPACK_IMPL(slaisnan)(SIN1,
                                 SIN2);
}

//-- slamch --------------------------------------------------------------------
FLOAT
LAPACK_DECL(slamch)(const char   *CMACH)
{
    DEBUG_LAPACK_STUB("slamch");
    return LAPACK_IMPL(slamch)(CMACH);
}

//-- slaswp --------------------------------------------------------------------
void
LAPACK_DECL(slaswp)(const INTEGER    *N,
                    FLOAT            *A,
                    const INTEGER    *LDA,
                    const INTEGER    *K1,
                    const INTEGER    *K2,
                    const INTEGER    *IPIV,
                    const INTEGER    *INCX)
{
    DEBUG_LAPACK_STUB("slaswp");
    LAPACK_IMPL(slaswp)(N,
                        A,
                        LDA,
                        K1,
                        K2,
                        IPIV,
                        INCX);
}

//-- spotf2 --------------------------------------------------------------------
void
LAPACK_DECL(spotf2)(const char       *UPLO,
                    const INTEGER    *N,
                    FLOAT            *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("spotf2");
    LAPACK_IMPL(spotf2)(UPLO,
                        N,
                        A,
                        LDA,
                        INFO);
}

//-- spotrf --------------------------------------------------------------------
void
LAPACK_DECL(spotrf)(const char       *UPLO,
                    const INTEGER    *N,
                    FLOAT            *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("spotrf");
    LAPACK_IMPL(spotrf)(UPLO,
                        N,
                        A,
                        LDA,
                        INFO);
}

//-- spotrs --------------------------------------------------------------------
void
LAPACK_DECL(spotrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const FLOAT      *A,
                    const INTEGER    *LDA,
                    FLOAT            *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("spotrs");
    LAPACK_IMPL(spotrs)(UPLO,
                        N,
                        NRHS,
                        A,
                        LDA,
                        B,
                        LDB,
                        INFO);
}

