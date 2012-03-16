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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zbbcsd");
    LAPACK_IMPL(zbbcsd)(JOBU1,
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
                        RWORK,
                        LRWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zbdsqr");
    LAPACK_IMPL(zbdsqr)(UPLO,
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
                        RWORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zcgesv");
    LAPACK_IMPL(zcgesv)(N,
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
                        RWORK,
                        ITER,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zcposv");
    LAPACK_IMPL(zcposv)(UPLO,
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
                        RWORK,
                        ITER,
                        INFO);
}

//-- zdrscl --------------------------------------------------------------------
void
LAPACK_DECL(zdrscl)(const INTEGER    *N,
                    const DOUBLE     *SA,
                    DOUBLE_COMPLEX   *SX,
                    const INTEGER    *INCX)
{
    DEBUG_LAPACK_STUB("zdrscl");
    LAPACK_IMPL(zdrscl)(N,
                        SA,
                        SX,
                        INCX);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgbbrd");
    LAPACK_IMPL(zgbbrd)(VECT,
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
                        RWORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgbcon");
    LAPACK_IMPL(zgbcon)(NORM,
                        N,
                        KL,
                        KU,
                        AB,
                        LDAB,
                        IPIV,
                        ANORM,
                        RCOND,
                        WORK,
                        RWORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgbequ");
    LAPACK_IMPL(zgbequ)(M,
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
                     INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgbequb");
    LAPACK_IMPL(zgbequb)(M,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgbrfs");
    LAPACK_IMPL(zgbrfs)(TRANS,
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
                        RWORK,
                        INFO);
}

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
                     INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgbrfsx");
    LAPACK_IMPL(zgbrfsx)(TRANS,
                         EQUED,
                         N,
                         KL,
                         KU,
                         NRHS,
                         AB,
                         LDAB,
                         AFB,
                         LDAFB,
                         IPIV,
                         R,
                         C,
                         B,
                         LDB,
                         X,
                         LDX,
                         RCOND,
                         BERR,
                         N_ERR_BNDS,
                         ERR_BNDS_NORM,
                         ERR_BNDS_COMP,
                         NPARAMS,
                         PARAMS,
                         WORK,
                         RWORK,
                         INFO);
}
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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zgbsv");
    LAPACK_IMPL(zgbsv)(N,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgbsvx");
    LAPACK_IMPL(zgbsvx)(FACT,
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
                        RWORK,
                        INFO);
}

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
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgbsvxx");
    LAPACK_IMPL(zgbsvxx)(FACT,
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
                         RPVGRW,
                         BERR,
                         N_ERR_BNDS,
                         ERR_BNDS_NORM,
                         ERR_BNDS_COMP,
                         NPARAMS,
                         PARAMS,
                         WORK,
                         RWORK,
                         INFO);
}
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgbtf2");
    LAPACK_IMPL(zgbtf2)(M,
                        N,
                        KL,
                        KU,
                        AB,
                        LDAB,
                        IPIV,
                        INFO);
}

//-- zgbtrf --------------------------------------------------------------------
void
LAPACK_DECL(zgbtrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *KL,
                    const INTEGER    *KU,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgbtrf");
    LAPACK_IMPL(zgbtrf)(M,
                        N,
                        KL,
                        KU,
                        AB,
                        LDAB,
                        IPIV,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgbtrs");
    LAPACK_IMPL(zgbtrs)(TRANS,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgebak");
    LAPACK_IMPL(zgebak)(JOB,
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

//-- zgebal --------------------------------------------------------------------
void
LAPACK_DECL(zgebal)(const char       *JOB,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *ILO,
                    INTEGER          *IHI,
                    DOUBLE           *SCALE,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgebal");
    LAPACK_IMPL(zgebal)(JOB,
                        N,
                        A,
                        LDA,
                        ILO,
                        IHI,
                        SCALE,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgebd2");
    LAPACK_IMPL(zgebd2)(M,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgebrd");
    LAPACK_IMPL(zgebrd)(M,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgecon");
    LAPACK_IMPL(zgecon)(NORM,
                        N,
                        A,
                        LDA,
                        ANORM,
                        RCOND,
                        WORK,
                        RWORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgeequ");
    LAPACK_IMPL(zgeequ)(M,
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
                     INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgeequb");
    LAPACK_IMPL(zgeequb)(M,
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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zgees");
    LAPACK_IMPL(zgees)(JOBVS,
                       SORT,
                       SELECT,
                       N,
                       A,
                       LDA,
                       SDIM,
                       W,
                       VS,
                       LDVS,
                       WORK,
                       LWORK,
                       RWORK,
                       BWORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgeesx");
    LAPACK_IMPL(zgeesx)(JOBVS,
                        SORT,
                        SELECT,
                        SENSE,
                        N,
                        A,
                        LDA,
                        SDIM,
                        W,
                        VS,
                        LDVS,
                        RCONDE,
                        RCONDV,
                        WORK,
                        LWORK,
                        RWORK,
                        BWORK,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zgeev");
    LAPACK_IMPL(zgeev)(JOBVL,
                       JOBVR,
                       N,
                       A,
                       LDA,
                       W,
                       VL,
                       LDVL,
                       VR,
                       LDVR,
                       WORK,
                       LWORK,
                       RWORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgeevx");
    LAPACK_IMPL(zgeevx)(BALANC,
                        JOBVL,
                        JOBVR,
                        SENSE,
                        N,
                        A,
                        LDA,
                        W,
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
                        RWORK,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zgegs");
    LAPACK_IMPL(zgegs)(JOBVSL,
                       JOBVSR,
                       N,
                       A,
                       LDA,
                       B,
                       LDB,
                       ALPHA,
                       BETA,
                       VSL,
                       LDVSL,
                       VSR,
                       LDVSR,
                       WORK,
                       LWORK,
                       RWORK,
                       INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zgegv");
    LAPACK_IMPL(zgegv)(JOBVL,
                       JOBVR,
                       N,
                       A,
                       LDA,
                       B,
                       LDB,
                       ALPHA,
                       BETA,
                       VL,
                       LDVL,
                       VR,
                       LDVR,
                       WORK,
                       LWORK,
                       RWORK,
                       INFO);
}

//-- zgehd2 --------------------------------------------------------------------
void
LAPACK_DECL(zgehd2)(const INTEGER    *N,
                    const INTEGER    *ILO,
                    const INTEGER    *IHI,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgehd2");
    LAPACK_IMPL(zgehd2)(N,
                        ILO,
                        IHI,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgehrd");
    LAPACK_IMPL(zgehrd)(N,
                        ILO,
                        IHI,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

//-- zgelq2 --------------------------------------------------------------------
void
LAPACK_DECL(zgelq2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgelq2");
    LAPACK_IMPL(zgelq2)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

//-- zgelqf --------------------------------------------------------------------
void
LAPACK_DECL(zgelqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgelqf");
    LAPACK_IMPL(zgelqf)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zgels");
    LAPACK_IMPL(zgels)(TRANS,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgelsd");
    LAPACK_IMPL(zgelsd)(M,
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
                        RWORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgelss");
    LAPACK_IMPL(zgelss)(M,
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
                        RWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgelsx");
    LAPACK_IMPL(zgelsx)(M,
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
                        RWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgelsy");
    LAPACK_IMPL(zgelsy)(M,
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
                        RWORK,
                        INFO);
}

//-- zgeql2 --------------------------------------------------------------------
void
LAPACK_DECL(zgeql2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgeql2");
    LAPACK_IMPL(zgeql2)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

//-- zgeqlf --------------------------------------------------------------------
void
LAPACK_DECL(zgeqlf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgeqlf");
    LAPACK_IMPL(zgeqlf)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgeqp3");
    LAPACK_IMPL(zgeqp3)(M,
                        N,
                        A,
                        LDA,
                        JPVT,
                        TAU,
                        WORK,
                        LWORK,
                        RWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgeqpf");
    LAPACK_IMPL(zgeqpf)(M,
                        N,
                        A,
                        LDA,
                        JPVT,
                        TAU,
                        WORK,
                        RWORK,
                        INFO);
}

//-- zgeqr2 --------------------------------------------------------------------
void
LAPACK_DECL(zgeqr2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgeqr2");
    LAPACK_IMPL(zgeqr2)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

//-- zgeqr2p -------------------------------------------------------------------
void
LAPACK_DECL(zgeqr2p)(const INTEGER    *M,
                     const INTEGER    *N,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     DOUBLE_COMPLEX   *TAU,
                     DOUBLE_COMPLEX   *WORK,
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgeqr2p");
    LAPACK_IMPL(zgeqr2p)(M,
                         N,
                         A,
                         LDA,
                         TAU,
                         WORK,
                         INFO);
}

//-- zgeqrf --------------------------------------------------------------------
void
LAPACK_DECL(zgeqrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgeqrf");
    LAPACK_IMPL(zgeqrf)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

//-- zgeqrfp -------------------------------------------------------------------
void
LAPACK_DECL(zgeqrfp)(const INTEGER    *M,
                     const INTEGER    *N,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     DOUBLE_COMPLEX   *TAU,
                     DOUBLE_COMPLEX   *WORK,
                     const INTEGER    *LWORK,
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgeqrfp");
    LAPACK_IMPL(zgeqrfp)(M,
                         N,
                         A,
                         LDA,
                         TAU,
                         WORK,
                         LWORK,
                         INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgerfs");
    LAPACK_IMPL(zgerfs)(TRANS,
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
                        RWORK,
                        INFO);
}

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
                     INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgerfsx");
    LAPACK_IMPL(zgerfsx)(TRANS,
                         EQUED,
                         N,
                         NRHS,
                         A,
                         LDA,
                         AF,
                         LDAF,
                         IPIV,
                         R,
                         C,
                         B,
                         LDB,
                         X,
                         LDX,
                         RCOND,
                         BERR,
                         N_ERR_BNDS,
                         ERR_BNDS_NORM,
                         ERR_BNDS_COMP,
                         NPARAMS,
                         PARAMS,
                         WORK,
                         RWORK,
                         INFO);
}
*/

//-- zgerq2 --------------------------------------------------------------------
void
LAPACK_DECL(zgerq2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgerq2");
    LAPACK_IMPL(zgerq2)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

//-- zgerqf --------------------------------------------------------------------
void
LAPACK_DECL(zgerqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgerqf");
    LAPACK_IMPL(zgerqf)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

//-- zgesc2 --------------------------------------------------------------------
void
LAPACK_DECL(zgesc2)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *RHS,
                    const INTEGER            *IPIV,
                    const INTEGER            *JPIV,
                    DOUBLE                   *SCALE)
{
    DEBUG_LAPACK_STUB("zgesc2");
    LAPACK_IMPL(zgesc2)(N,
                        A,
                        LDA,
                        RHS,
                        IPIV,
                        JPIV,
                        SCALE);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgesdd");
    LAPACK_IMPL(zgesdd)(JOBZ,
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
                        RWORK,
                        IWORK,
                        INFO);
}

//-- zgesv ---------------------------------------------------------------------
void
LAPACK_DECL(zgesv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   INTEGER              *IPIV,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zgesv");
    LAPACK_IMPL(zgesv)(N,
                       NRHS,
                       A,
                       LDA,
                       IPIV,
                       B,
                       LDB,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgesvd");
    LAPACK_IMPL(zgesvd)(JOBU,
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
                        RWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgesvx");
    LAPACK_IMPL(zgesvx)(FACT,
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
                        RWORK,
                        INFO);
}

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
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgesvxx");
    LAPACK_IMPL(zgesvxx)(FACT,
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
                         RPVGRW,
                         BERR,
                         N_ERR_BNDS,
                         ERR_BNDS_NORM,
                         ERR_BNDS_COMP,
                         NPARAMS,
                         PARAMS,
                         WORK,
                         RWORK,
                         INFO);
}
*/
//-- zgetc2 --------------------------------------------------------------------
void
LAPACK_DECL(zgetc2)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *JPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgetc2");
    LAPACK_IMPL(zgetc2)(N,
                        A,
                        LDA,
                        IPIV,
                        JPIV,
                        INFO);
}

//-- zgetf2 --------------------------------------------------------------------
void
LAPACK_DECL(zgetf2)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgetf2");
    LAPACK_IMPL(zgetf2)(M,
                        N,
                        A,
                        LDA,
                        IPIV,
                        INFO);
}

//-- zgetrf --------------------------------------------------------------------
/*
void
LAPACK_DECL(zgetrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgetrf");
    LAPACK_IMPL(zgetrf)(M,
                        N,
                        A,
                        LDA,
                        IPIV,
                        INFO);
}
*/

//-- zgetri --------------------------------------------------------------------
void
LAPACK_DECL(zgetri)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgetri");
    LAPACK_IMPL(zgetri)(N,
                        A,
                        LDA,
                        IPIV,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgetrs");
    LAPACK_IMPL(zgetrs)(TRANS,
                        N,
                        NRHS,
                        A,
                        LDA,
                        IPIV,
                        B,
                        LDB,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zggbak");
    LAPACK_IMPL(zggbak)(JOB,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zggbal");
    LAPACK_IMPL(zggbal)(JOB,
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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zgges");
    LAPACK_IMPL(zgges)(JOBVSL,
                       JOBVSR,
                       SORT,
                       SELCTG,
                       N,
                       A,
                       LDA,
                       B,
                       LDB,
                       SDIM,
                       ALPHA,
                       BETA,
                       VSL,
                       LDVSL,
                       VSR,
                       LDVSR,
                       WORK,
                       LWORK,
                       RWORK,
                       BWORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zggesx");
    LAPACK_IMPL(zggesx)(JOBVSL,
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
                        ALPHA,
                        BETA,
                        VSL,
                        LDVSL,
                        VSR,
                        LDVSR,
                        RCONDE,
                        RCONDV,
                        WORK,
                        LWORK,
                        RWORK,
                        IWORK,
                        LIWORK,
                        BWORK,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zggev");
    LAPACK_IMPL(zggev)(JOBVL,
                       JOBVR,
                       N,
                       A,
                       LDA,
                       B,
                       LDB,
                       ALPHA,
                       BETA,
                       VL,
                       LDVL,
                       VR,
                       LDVR,
                       WORK,
                       LWORK,
                       RWORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zggevx");
    LAPACK_IMPL(zggevx)(BALANC,
                        JOBVL,
                        JOBVR,
                        SENSE,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        ALPHA,
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
                        RWORK,
                        IWORK,
                        BWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zggglm");
    LAPACK_IMPL(zggglm)(N,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgghrd");
    LAPACK_IMPL(zgghrd)(COMPQ,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgglse");
    LAPACK_IMPL(zgglse)(M,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zggqrf");
    LAPACK_IMPL(zggqrf)(N,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zggrqf");
    LAPACK_IMPL(zggrqf)(M,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zggsvd");
    LAPACK_IMPL(zggsvd)(JOBU,
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
                        RWORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zggsvp");
    LAPACK_IMPL(zggsvp)(JOBU,
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
                        RWORK,
                        TAU,
                        WORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgtcon");
    LAPACK_IMPL(zgtcon)(NORM,
                        N,
                        DL,
                        D,
                        DU,
                        DU2,
                        IPIV,
                        ANORM,
                        RCOND,
                        WORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgtrfs");
    LAPACK_IMPL(zgtrfs)(TRANS,
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
                        RWORK,
                        INFO);
}

//-- zgtsv ---------------------------------------------------------------------
void
LAPACK_DECL(zgtsv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *DL,
                   DOUBLE_COMPLEX       *D,
                   DOUBLE_COMPLEX       *DU,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zgtsv");
    LAPACK_IMPL(zgtsv)(N,
                       NRHS,
                       DL,
                       D,
                       DU,
                       B,
                       LDB,
                       INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgtsvx");
    LAPACK_IMPL(zgtsvx)(FACT,
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
                        RWORK,
                        INFO);
}

//-- zgttrf --------------------------------------------------------------------
void
LAPACK_DECL(zgttrf)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *DL,
                    DOUBLE_COMPLEX   *D,
                    DOUBLE_COMPLEX   *DU,
                    DOUBLE_COMPLEX   *DU2,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zgttrf");
    LAPACK_IMPL(zgttrf)(N,
                        DL,
                        D,
                        DU,
                        DU2,
                        IPIV,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zgttrs");
    LAPACK_IMPL(zgttrs)(TRANS,
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
                    const INTEGER            *LDB)
{
    DEBUG_LAPACK_STUB("zgtts2");
    LAPACK_IMPL(zgtts2)(ITRANS,
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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zhbev");
    LAPACK_IMPL(zhbev)(JOBZ,
                       UPLO,
                       N,
                       KD,
                       AB,
                       LDAB,
                       W,
                       Z,
                       LDZ,
                       WORK,
                       RWORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhbevd");
    LAPACK_IMPL(zhbevd)(JOBZ,
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
                        RWORK,
                        LRWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhbevx");
    LAPACK_IMPL(zhbevx)(JOBZ,
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
                        RWORK,
                        IWORK,
                        IFAIL,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zhbgst");
    LAPACK_IMPL(zhbgst)(VECT,
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
                        RWORK,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zhbgv");
    LAPACK_IMPL(zhbgv)(JOBZ,
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
                       RWORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhbgvd");
    LAPACK_IMPL(zhbgvd)(JOBZ,
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
                        RWORK,
                        LRWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhbgvx");
    LAPACK_IMPL(zhbgvx)(JOBZ,
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
                        RWORK,
                        IWORK,
                        IFAIL,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhbtrd");
    LAPACK_IMPL(zhbtrd)(VECT,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zhecon");
    LAPACK_IMPL(zhecon)(UPLO,
                        N,
                        A,
                        LDA,
                        IPIV,
                        ANORM,
                        RCOND,
                        WORK,
                        INFO);
}

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
                     INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zheequb");
    LAPACK_IMPL(zheequb)(UPLO,
                         N,
                         A,
                         LDA,
                         S,
                         SCOND,
                         AMAX,
                         WORK,
                         INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zheev");
    LAPACK_IMPL(zheev)(JOBZ,
                       UPLO,
                       N,
                       A,
                       LDA,
                       W,
                       WORK,
                       LWORK,
                       RWORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zheevd");
    LAPACK_IMPL(zheevd)(JOBZ,
                        UPLO,
                        N,
                        A,
                        LDA,
                        W,
                        WORK,
                        LWORK,
                        RWORK,
                        LRWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zheevr");
    LAPACK_IMPL(zheevr)(JOBZ,
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
                        RWORK,
                        LRWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zheevx");
    LAPACK_IMPL(zheevx)(JOBZ,
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
                        RWORK,
                        IWORK,
                        IFAIL,
                        INFO);
}

//-- zhegs2 --------------------------------------------------------------------
void
LAPACK_DECL(zhegs2)(const INTEGER            *ITYPE,
                    const char               *UPLO,
                    const INTEGER            *N,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zhegs2");
    LAPACK_IMPL(zhegs2)(ITYPE,
                        UPLO,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        INFO);
}

//-- zhegst --------------------------------------------------------------------
void
LAPACK_DECL(zhegst)(const INTEGER            *ITYPE,
                    const char               *UPLO,
                    const INTEGER            *N,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zhegst");
    LAPACK_IMPL(zhegst)(ITYPE,
                        UPLO,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zhegv");
    LAPACK_IMPL(zhegv)(ITYPE,
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
                       RWORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhegvd");
    LAPACK_IMPL(zhegvd)(ITYPE,
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
                        RWORK,
                        LRWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhegvx");
    LAPACK_IMPL(zhegvx)(ITYPE,
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
                        RWORK,
                        IWORK,
                        IFAIL,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zherfs");
    LAPACK_IMPL(zherfs)(UPLO,
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
                        RWORK,
                        INFO);
}

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
                     INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zherfsx");
    LAPACK_IMPL(zherfsx)(UPLO,
                         EQUED,
                         N,
                         NRHS,
                         A,
                         LDA,
                         AF,
                         LDAF,
                         IPIV,
                         S,
                         B,
                         LDB,
                         X,
                         LDX,
                         RCOND,
                         BERR,
                         N_ERR_BNDS,
                         ERR_BNDS_NORM,
                         ERR_BNDS_COMP,
                         NPARAMS,
                         PARAMS,
                         WORK,
                         RWORK,
                         INFO);
}
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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zhesv");
    LAPACK_IMPL(zhesv)(UPLO,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zhesvx");
    LAPACK_IMPL(zhesvx)(FACT,
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
                        RWORK,
                        INFO);
}

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
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhesvxx");
    LAPACK_IMPL(zhesvxx)(FACT,
                         UPLO,
                         N,
                         NRHS,
                         A,
                         LDA,
                         AF,
                         LDAF,
                         IPIV,
                         EQUED,
                         S,
                         B,
                         LDB,
                         X,
                         LDX,
                         RCOND,
                         RPVGRW,
                         BERR,
                         N_ERR_BNDS,
                         ERR_BNDS_NORM,
                         ERR_BNDS_COMP,
                         NPARAMS,
                         PARAMS,
                         WORK,
                         RWORK,
                         INFO);
}
*/
//-- zheswapr ------------------------------------------------------------------
void
LAPACK_DECL(zheswapr)(const char           *UPLO,
                      const INTEGER        *N,
                      DOUBLE_COMPLEX       *A,
                      const INTEGER        *LDA,
                      const INTEGER        *I1,
                      const INTEGER        *I2)
{
    DEBUG_LAPACK_STUB("zheswapr");
    LAPACK_IMPL(zheswapr)(UPLO,
                          N,
                          A,
                          LDA,
                          I1,
                          I2);
}

//-- zhetd2 --------------------------------------------------------------------
void
LAPACK_DECL(zhetd2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *TAU,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhetd2");
    LAPACK_IMPL(zhetd2)(UPLO,
                        N,
                        A,
                        LDA,
                        D,
                        E,
                        TAU,
                        INFO);
}

//-- zhetf2 --------------------------------------------------------------------
void
LAPACK_DECL(zhetf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhetf2");
    LAPACK_IMPL(zhetf2)(UPLO,
                        N,
                        A,
                        LDA,
                        IPIV,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhetrd");
    LAPACK_IMPL(zhetrd)(UPLO,
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

//-- zhetrf --------------------------------------------------------------------
void
LAPACK_DECL(zhetrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhetrf");
    LAPACK_IMPL(zhetrf)(UPLO,
                        N,
                        A,
                        LDA,
                        IPIV,
                        WORK,
                        LWORK,
                        INFO);
}

//-- zhetri --------------------------------------------------------------------
void
LAPACK_DECL(zhetri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhetri");
    LAPACK_IMPL(zhetri)(UPLO,
                        N,
                        A,
                        LDA,
                        IPIV,
                        WORK,
                        INFO);
}

//-- zhetri2 -------------------------------------------------------------------
void
LAPACK_DECL(zhetri2)(const char       *UPLO,
                     const INTEGER    *N,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE_COMPLEX   *WORK,
                     const INTEGER    *LWORK,
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhetri2");
    LAPACK_IMPL(zhetri2)(UPLO,
                         N,
                         A,
                         LDA,
                         IPIV,
                         WORK,
                         LWORK,
                         INFO);
}

//-- zhetri2x ------------------------------------------------------------------
void
LAPACK_DECL(zhetri2x)(const char           *UPLO,
                      const INTEGER        *N,
                      DOUBLE_COMPLEX       *A,
                      const INTEGER        *LDA,
                      const INTEGER        *IPIV,
                      DOUBLE_COMPLEX       *WORK,
                      const INTEGER        *NB,
                      INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zhetri2x");
    LAPACK_IMPL(zhetri2x)(UPLO,
                          N,
                          A,
                          LDA,
                          IPIV,
                          WORK,
                          NB,
                          INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zhetrs");
    LAPACK_IMPL(zhetrs)(UPLO,
                        N,
                        NRHS,
                        A,
                        LDA,
                        IPIV,
                        B,
                        LDB,
                        INFO);
}

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
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhetrs2");
    LAPACK_IMPL(zhetrs2)(UPLO,
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
                   DOUBLE_COMPLEX           *C)
{
    DEBUG_LAPACK_STUB("zhfrk");
    LAPACK_IMPL(zhfrk)(TRANSR,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhgeqz");
    LAPACK_IMPL(zhgeqz)(JOB,
                        COMPQ,
                        COMPZ,
                        N,
                        ILO,
                        IHI,
                        H,
                        LDH,
                        T,
                        LDT,
                        ALPHA,
                        BETA,
                        Q,
                        LDQ,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        RWORK,
                        INFO);
}

//-- zhpcon --------------------------------------------------------------------
void
LAPACK_DECL(zhpcon)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    const INTEGER            *IPIV,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zhpcon");
    LAPACK_IMPL(zhpcon)(UPLO,
                        N,
                        AP,
                        IPIV,
                        ANORM,
                        RCOND,
                        WORK,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zhpev");
    LAPACK_IMPL(zhpev)(JOBZ,
                       UPLO,
                       N,
                       AP,
                       W,
                       Z,
                       LDZ,
                       WORK,
                       RWORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhpevd");
    LAPACK_IMPL(zhpevd)(JOBZ,
                        UPLO,
                        N,
                        AP,
                        W,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        RWORK,
                        LRWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhpevx");
    LAPACK_IMPL(zhpevx)(JOBZ,
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
                        RWORK,
                        IWORK,
                        IFAIL,
                        INFO);
}

//-- zhpgst --------------------------------------------------------------------
void
LAPACK_DECL(zhpgst)(const INTEGER            *ITYPE,
                    const char               *UPLO,
                    const INTEGER            *N,
                    DOUBLE_COMPLEX           *AP,
                    const DOUBLE_COMPLEX     *BP,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zhpgst");
    LAPACK_IMPL(zhpgst)(ITYPE,
                        UPLO,
                        N,
                        AP,
                        BP,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zhpgv");
    LAPACK_IMPL(zhpgv)(ITYPE,
                       JOBZ,
                       UPLO,
                       N,
                       AP,
                       BP,
                       W,
                       Z,
                       LDZ,
                       WORK,
                       RWORK,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhpgvd");
    LAPACK_IMPL(zhpgvd)(ITYPE,
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
                        RWORK,
                        LRWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhpgvx");
    LAPACK_IMPL(zhpgvx)(ITYPE,
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
                        RWORK,
                        IWORK,
                        IFAIL,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zhprfs");
    LAPACK_IMPL(zhprfs)(UPLO,
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
                        RWORK,
                        INFO);
}

//-- zhpsv ---------------------------------------------------------------------
void
LAPACK_DECL(zhpsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *AP,
                   INTEGER              *IPIV,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zhpsv");
    LAPACK_IMPL(zhpsv)(UPLO,
                       N,
                       NRHS,
                       AP,
                       IPIV,
                       B,
                       LDB,
                       INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zhpsvx");
    LAPACK_IMPL(zhpsvx)(FACT,
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
                        RWORK,
                        INFO);
}

//-- zhptrd --------------------------------------------------------------------
void
LAPACK_DECL(zhptrd)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *TAU,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhptrd");
    LAPACK_IMPL(zhptrd)(UPLO,
                        N,
                        AP,
                        D,
                        E,
                        TAU,
                        INFO);
}

//-- zhptrf --------------------------------------------------------------------
void
LAPACK_DECL(zhptrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhptrf");
    LAPACK_IMPL(zhptrf)(UPLO,
                        N,
                        AP,
                        IPIV,
                        INFO);
}

//-- zhptri --------------------------------------------------------------------
void
LAPACK_DECL(zhptri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    const INTEGER    *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhptri");
    LAPACK_IMPL(zhptri)(UPLO,
                        N,
                        AP,
                        IPIV,
                        WORK,
                        INFO);
}

//-- zhptrs --------------------------------------------------------------------
void
LAPACK_DECL(zhptrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    const INTEGER            *IPIV,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zhptrs");
    LAPACK_IMPL(zhptrs)(UPLO,
                        N,
                        NRHS,
                        AP,
                        IPIV,
                        B,
                        LDB,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zhsein");
    LAPACK_IMPL(zhsein)(SIDE,
                        EIGSRC,
                        INITV,
                        SELECT,
                        N,
                        H,
                        LDH,
                        W,
                        VL,
                        LDVL,
                        VR,
                        LDVR,
                        MM,
                        M,
                        WORK,
                        RWORK,
                        IFAILL,
                        IFAILR,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zhseqr");
    LAPACK_IMPL(zhseqr)(JOB,
                        COMPZ,
                        N,
                        ILO,
                        IHI,
                        H,
                        LDH,
                        W,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        INFO);
}

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
                       const INTEGER            *INCY)
{
    DEBUG_LAPACK_STUB("zla_gbamv");
    LAPACK_IMPL(zla_gbamv)(TRANS,
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
                           const DOUBLE             *RWORK)
{
    DEBUG_LAPACK_STUB("zla_gbrcond_c");
    return LAPACK_IMPL(zla_gbrcond_c)(TRANS,
                                      N,
                                      KL,
                                      KU,
                                      AB,
                                      LDAB,
                                      AFB,
                                      LDAFB,
                                      IPIV,
                                      C,
                                      CAPPLY,
                                      INFO,
                                      WORK,
                                      RWORK);
}

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
                           const DOUBLE             *RWORK)
{
    DEBUG_LAPACK_STUB("zla_gbrcond_x");
    return LAPACK_IMPL(zla_gbrcond_x)(TRANS,
                                      N,
                                      KL,
                                      KU,
                                      AB,
                                      LDAB,
                                      AFB,
                                      LDAFB,
                                      IPIV,
                                      X,
                                      INFO,
                                      WORK,
                                      RWORK);
}

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
                                 INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zla_gbrfsx_extended");
    LAPACK_IMPL(zla_gbrfsx_extended)(PREC_TYPE,
                                     TRANS_TYPE,
                                     N,
                                     KL,
                                     KU,
                                     NRHS,
                                     AB,
                                     LDAB,
                                     AFB,
                                     LDAFB,
                                     IPIV,
                                     COLEQU,
                                     C,
                                     B,
                                     LDB,
                                     Y,
                                     LDY,
                                     BERR_OUT,
                                     N_NORMS,
                                     ERR_BNDS_NORM,
                                     ERR_BNDS_COMP,
                                     RES,
                                     AYB,
                                     DY,
                                     Y_TAIL,
                                     RCOND,
                                     ITHRESH,
                                     RTHRESH,
                                     DZ_UB,
                                     IGNORE_CWISE,
                                     INFO);
}
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
                          const INTEGER            *LDAFB)
{
    DEBUG_LAPACK_STUB("zla_gbrpvgrw");
    return LAPACK_IMPL(zla_gbrpvgrw)(N,
                                     KL,
                                     KU,
                                     NCOLS,
                                     AB,
                                     LDAB,
                                     AFB,
                                     LDAFB);
}

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
                       const INTEGER            *INCY)
{
    DEBUG_LAPACK_STUB("zla_geamv");
    LAPACK_IMPL(zla_geamv)(TRANS,
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
                           const DOUBLE             *RWORK)
{
    DEBUG_LAPACK_STUB("zla_gercond_c");
    return LAPACK_IMPL(zla_gercond_c)(TRANS,
                                      N,
                                      A,
                                      LDA,
                                      AF,
                                      LDAF,
                                      IPIV,
                                      C,
                                      CAPPLY,
                                      INFO,
                                      WORK,
                                      RWORK);
}

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
                           const DOUBLE             *RWORK)
{
    DEBUG_LAPACK_STUB("zla_gercond_x");
    return LAPACK_IMPL(zla_gercond_x)(TRANS,
                                      N,
                                      A,
                                      LDA,
                                      AF,
                                      LDAF,
                                      IPIV,
                                      X,
                                      INFO,
                                      WORK,
                                      RWORK);
}

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
                                 INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zla_gerfsx_extended");
    LAPACK_IMPL(zla_gerfsx_extended)(PREC_TYPE,
                                     TRANS_TYPE,
                                     N,
                                     NRHS,
                                     A,
                                     LDA,
                                     AF,
                                     LDAF,
                                     IPIV,
                                     COLEQU,
                                     C,
                                     B,
                                     LDB,
                                     Y,
                                     LDY,
                                     BERR_OUT,
                                     N_NORMS,
                                     ERRS_N,
                                     ERRS_C,
                                     RES,
                                     AYB,
                                     DY,
                                     Y_TAIL,
                                     RCOND,
                                     ITHRESH,
                                     RTHRESH,
                                     DZ_UB,
                                     IGNORE_CWISE,
                                     INFO);
}
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
                       const INTEGER            *INCY)
{
    DEBUG_LAPACK_STUB("zla_heamv");
    LAPACK_IMPL(zla_heamv)(UPLO,
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
                           const DOUBLE             *RWORK)
{
    DEBUG_LAPACK_STUB("zla_hercond_c");
    return LAPACK_IMPL(zla_hercond_c)(UPLO,
                                      N,
                                      A,
                                      LDA,
                                      AF,
                                      LDAF,
                                      IPIV,
                                      C,
                                      CAPPLY,
                                      INFO,
                                      WORK,
                                      RWORK);
}

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
                           const DOUBLE             *RWORK)
{
    DEBUG_LAPACK_STUB("zla_hercond_x");
    return LAPACK_IMPL(zla_hercond_x)(UPLO,
                                      N,
                                      A,
                                      LDA,
                                      AF,
                                      LDAF,
                                      IPIV,
                                      X,
                                      INFO,
                                      WORK,
                                      RWORK);
}

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
                                 INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zla_herfsx_extended");
    LAPACK_IMPL(zla_herfsx_extended)(PREC_TYPE,
                                     UPLO,
                                     N,
                                     NRHS,
                                     A,
                                     LDA,
                                     AF,
                                     LDAF,
                                     IPIV,
                                     COLEQU,
                                     C,
                                     B,
                                     LDB,
                                     Y,
                                     LDY,
                                     BERR_OUT,
                                     N_NORMS,
                                     ERR_BNDS_NORM,
                                     ERR_BNDS_COMP,
                                     RES,
                                     AYB,
                                     DY,
                                     Y_TAIL,
                                     RCOND,
                                     ITHRESH,
                                     RTHRESH,
                                     DZ_UB,
                                     IGNORE_CWISE,
                                     INFO);
}
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
                          const DOUBLE             *WORK)
{
    DEBUG_LAPACK_STUB("zla_herpvgrw");
    return LAPACK_IMPL(zla_herpvgrw)(UPLO,
                                     N,
                                     INFO,
                                     A,
                                     LDA,
                                     AF,
                                     LDAF,
                                     IPIV,
                                     WORK);
}

//-- zla_lin_berr --------------------------------------------------------------
void
LAPACK_DECL(zla_lin_berr)(const INTEGER            *N,
                          const INTEGER            *NZ,
                          const INTEGER            *NRHS,
                          const DOUBLE_COMPLEX     *RES,
                          const DOUBLE             *AYB,
                          DOUBLE                   *BERR)
{
    DEBUG_LAPACK_STUB("zla_lin_berr");
    LAPACK_IMPL(zla_lin_berr)(N,
                              NZ,
                              NRHS,
                              RES,
                              AYB,
                              BERR);
}

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
                           const DOUBLE             *RWORK)
{
    DEBUG_LAPACK_STUB("zla_porcond_c");
    return LAPACK_IMPL(zla_porcond_c)(UPLO,
                                      N,
                                      A,
                                      LDA,
                                      AF,
                                      LDAF,
                                      C,
                                      CAPPLY,
                                      INFO,
                                      WORK,
                                      RWORK);
}

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
                           const DOUBLE             *RWORK)
{
    DEBUG_LAPACK_STUB("zla_porcond_x");
    return LAPACK_IMPL(zla_porcond_x)(UPLO,
                                      N,
                                      A,
                                      LDA,
                                      AF,
                                      LDAF,
                                      X,
                                      INFO,
                                      WORK,
                                      RWORK);
}

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
                                 INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zla_porfsx_extended");
    LAPACK_IMPL(zla_porfsx_extended)(PREC_TYPE,
                                     UPLO,
                                     N,
                                     NRHS,
                                     A,
                                     LDA,
                                     AF,
                                     LDAF,
                                     COLEQU,
                                     C,
                                     B,
                                     LDB,
                                     Y,
                                     LDY,
                                     BERR_OUT,
                                     N_NORMS,
                                     ERR_BNDS_NORM,
                                     ERR_BNDS_COMP,
                                     RES,
                                     AYB,
                                     DY,
                                     Y_TAIL,
                                     RCOND,
                                     ITHRESH,
                                     RTHRESH,
                                     DZ_UB,
                                     IGNORE_CWISE,
                                     INFO);
}
*/

//-- zla_porpvgrw --------------------------------------------------------------
DOUBLE
LAPACK_DECL(zla_porpvgrw)(const char               *UPLO,
                          const INTEGER            *NCOLS,
                          const DOUBLE_COMPLEX     *A,
                          const INTEGER            *LDA,
                          const DOUBLE_COMPLEX     *AF,
                          const INTEGER            *LDAF,
                          const DOUBLE             *WORK)
{
    DEBUG_LAPACK_STUB("zla_porpvgrw");
    return LAPACK_IMPL(zla_porpvgrw)(UPLO,
                                     NCOLS,
                                     A,
                                     LDA,
                                     AF,
                                     LDAF,
                                     WORK);
}

//-- zla_rpvgrw ----------------------------------------------------------------
DOUBLE
LAPACK_DECL(zla_rpvgrw)(const INTEGER            *N,
                        const INTEGER            *NCOLS,
                        const DOUBLE_COMPLEX     *A,
                        const INTEGER            *LDA,
                        const DOUBLE_COMPLEX     *AF,
                        const INTEGER            *LDAF)
{
    DEBUG_LAPACK_STUB("zla_rpvgrw");
    return LAPACK_IMPL(zla_rpvgrw)(N,
                                   NCOLS,
                                   A,
                                   LDA,
                                   AF,
                                   LDAF);
}

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
                       const INTEGER            *INCY)
{
    DEBUG_LAPACK_STUB("zla_syamv");
    LAPACK_IMPL(zla_syamv)(UPLO,
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
                           const DOUBLE             *RWORK)
{
    DEBUG_LAPACK_STUB("zla_syrcond_c");
    return LAPACK_IMPL(zla_syrcond_c)(UPLO,
                                      N,
                                      A,
                                      LDA,
                                      AF,
                                      LDAF,
                                      IPIV,
                                      C,
                                      CAPPLY,
                                      INFO,
                                      WORK,
                                      RWORK);
}

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
                           const DOUBLE             *RWORK)
{
    DEBUG_LAPACK_STUB("zla_syrcond_x");
    return LAPACK_IMPL(zla_syrcond_x)(UPLO,
                                      N,
                                      A,
                                      LDA,
                                      AF,
                                      LDAF,
                                      IPIV,
                                      X,
                                      INFO,
                                      WORK,
                                      RWORK);
}

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
                                 INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zla_syrfsx_extended");
    LAPACK_IMPL(zla_syrfsx_extended)(PREC_TYPE,
                                     UPLO,
                                     N,
                                     NRHS,
                                     A,
                                     LDA,
                                     AF,
                                     LDAF,
                                     IPIV,
                                     COLEQU,
                                     C,
                                     B,
                                     LDB,
                                     Y,
                                     LDY,
                                     BERR_OUT,
                                     N_NORMS,
                                     ERR_BNDS_NORM,
                                     ERR_BNDS_COMP,
                                     RES,
                                     AYB,
                                     DY,
                                     Y_TAIL,
                                     RCOND,
                                     ITHRESH,
                                     RTHRESH,
                                     DZ_UB,
                                     IGNORE_CWISE,
                                     INFO);
}
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
                          const DOUBLE             *WORK)
{
    DEBUG_LAPACK_STUB("zla_syrpvgrw");
    return LAPACK_IMPL(zla_syrpvgrw)(UPLO,
                                     N,
                                     INFO,
                                     A,
                                     LDA,
                                     AF,
                                     LDAF,
                                     IPIV,
                                     WORK);
}

//-- zla_wwaddw ----------------------------------------------------------------
void
LAPACK_DECL(zla_wwaddw)(const INTEGER            *N,
                        DOUBLE_COMPLEX           *X,
                        DOUBLE_COMPLEX           *Y,
                        const DOUBLE_COMPLEX     *W)
{
    DEBUG_LAPACK_STUB("zla_wwaddw");
    LAPACK_IMPL(zla_wwaddw)(N,
                            X,
                            Y,
                            W);
}

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
                    const INTEGER    *LDY)
{
    DEBUG_LAPACK_STUB("zlabrd");
    LAPACK_IMPL(zlabrd)(M,
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

//-- zlacgv --------------------------------------------------------------------
void
LAPACK_DECL(zlacgv)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *INCX)
{
    DEBUG_LAPACK_STUB("zlacgv");
    LAPACK_IMPL(zlacgv)(N,
                        X,
                        INCX);
}

//-- zlacn2 --------------------------------------------------------------------
void
LAPACK_DECL(zlacn2)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *V,
                    DOUBLE_COMPLEX   *X,
                    DOUBLE           *EST,
                    INTEGER          *KASE,
                    INTEGER          *ISAVE)
{
    DEBUG_LAPACK_STUB("zlacn2");
    LAPACK_IMPL(zlacn2)(N,
                        V,
                        X,
                        EST,
                        KASE,
                        ISAVE);
}

//-- zlacon --------------------------------------------------------------------
void
LAPACK_DECL(zlacon)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *V,
                    DOUBLE_COMPLEX   *X,
                    DOUBLE           *EST,
                    INTEGER          *KASE)
{
    DEBUG_LAPACK_STUB("zlacon");
    LAPACK_IMPL(zlacon)(N,
                        V,
                        X,
                        EST,
                        KASE);
}

//-- zlacp2 --------------------------------------------------------------------
void
LAPACK_DECL(zlacp2)(const char       *UPLO,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *B,
                    const INTEGER    *LDB)
{
    DEBUG_LAPACK_STUB("zlacp2");
    LAPACK_IMPL(zlacp2)(UPLO,
                        M,
                        N,
                        A,
                        LDA,
                        B,
                        LDB);
}

//-- zlacpy --------------------------------------------------------------------
void
LAPACK_DECL(zlacpy)(const char               *UPLO,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB)
{
    DEBUG_LAPACK_STUB("zlacpy");
    LAPACK_IMPL(zlacpy)(UPLO,
                        M,
                        N,
                        A,
                        LDA,
                        B,
                        LDB);
}

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
                    DOUBLE                   *RWORK)
{
    DEBUG_LAPACK_STUB("zlacrm");
    LAPACK_IMPL(zlacrm)(M,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        C,
                        LDC,
                        RWORK);
}

//-- zlacrt --------------------------------------------------------------------
void
LAPACK_DECL(zlacrt)(const INTEGER            *N,
                    DOUBLE_COMPLEX           *CX,
                    const INTEGER            *INCX,
                    DOUBLE_COMPLEX           *CY,
                    const INTEGER            *INCY,
                    const DOUBLE_COMPLEX     *C,
                    const DOUBLE_COMPLEX     *S)
{
    DEBUG_LAPACK_STUB("zlacrt");
    LAPACK_IMPL(zlacrt)(N,
                        CX,
                        INCX,
                        CY,
                        INCY,
                        C,
                        S);
}

//-- zladiv --------------------------------------------------------------------
/*
UNKNOWN
LAPACK_DECL(zladiv)(const DOUBLE_COMPLEX     *X,
                    const DOUBLE_COMPLEX     *Y)
{
    DEBUG_LAPACK_STUB("zladiv");
    return LAPACK_IMPL(zladiv)(X,
                               Y);
}
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zlaed0");
    LAPACK_IMPL(zlaed0)(QSIZ,
                        N,
                        D,
                        E,
                        Q,
                        LDQ,
                        QSTORE,
                        LDQS,
                        RWORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zlaed7");
    LAPACK_IMPL(zlaed7)(N,
                        CUTPNT,
                        QSIZ,
                        TLVLS,
                        CURLVL,
                        CURPBM,
                        D,
                        Q,
                        LDQ,
                        RHO,
                        INDXQ,
                        QSTORE,
                        QPTR,
                        PRMPTR,
                        PERM,
                        GIVPTR,
                        GIVCOL,
                        GIVNUM,
                        WORK,
                        RWORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zlaed8");
    LAPACK_IMPL(zlaed8)(K,
                        N,
                        QSIZ,
                        Q,
                        LDQ,
                        D,
                        RHO,
                        CUTPNT,
                        Z,
                        DLAMDA,
                        Q2,
                        LDQ2,
                        W,
                        INDXP,
                        INDX,
                        INDXQ,
                        PERM,
                        GIVPTR,
                        GIVCOL,
                        GIVNUM,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zlaein");
    LAPACK_IMPL(zlaein)(RIGHTV,
                        NOINIT,
                        N,
                        H,
                        LDH,
                        W,
                        V,
                        B,
                        LDB,
                        RWORK,
                        EPS3,
                        SMLNUM,
                        INFO);
}

//-- zlaesy --------------------------------------------------------------------
void
LAPACK_DECL(zlaesy)(const DOUBLE_COMPLEX     *A,
                    const DOUBLE_COMPLEX     *B,
                    const DOUBLE_COMPLEX     *C,
                    DOUBLE_COMPLEX           *RT1,
                    DOUBLE_COMPLEX           *RT2,
                    DOUBLE_COMPLEX           *EVSCAL,
                    DOUBLE_COMPLEX           *CS1,
                    DOUBLE_COMPLEX           *SN1)
{
    DEBUG_LAPACK_STUB("zlaesy");
    LAPACK_IMPL(zlaesy)(A,
                        B,
                        C,
                        RT1,
                        RT2,
                        EVSCAL,
                        CS1,
                        SN1);
}

//-- zlaev2 --------------------------------------------------------------------
void
LAPACK_DECL(zlaev2)(const DOUBLE_COMPLEX     *A,
                    const DOUBLE_COMPLEX     *B,
                    const DOUBLE_COMPLEX     *C,
                    DOUBLE                   *RT1,
                    DOUBLE                   *RT2,
                    DOUBLE                   *CS1,
                    DOUBLE_COMPLEX           *SN1)
{
    DEBUG_LAPACK_STUB("zlaev2");
    LAPACK_IMPL(zlaev2)(A,
                        B,
                        C,
                        RT1,
                        RT2,
                        CS1,
                        SN1);
}

//-- zlag2c --------------------------------------------------------------------
void
LAPACK_DECL(zlag2c)(const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    FLOAT_COMPLEX            *SA,
                    const INTEGER            *LDSA,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zlag2c");
    LAPACK_IMPL(zlag2c)(M,
                        N,
                        A,
                        LDA,
                        SA,
                        LDSA,
                        INFO);
}

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
                    DOUBLE_COMPLEX           *SNQ)
{
    DEBUG_LAPACK_STUB("zlags2");
    LAPACK_IMPL(zlags2)(UPPER,
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
                    const INTEGER            *LDB)
{
    DEBUG_LAPACK_STUB("zlagtm");
    LAPACK_IMPL(zlagtm)(TRANS,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zlahef");
    LAPACK_IMPL(zlahef)(UPLO,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zlahqr");
    LAPACK_IMPL(zlahqr)(WANTT,
                        WANTZ,
                        N,
                        ILO,
                        IHI,
                        H,
                        LDH,
                        W,
                        ILOZ,
                        IHIZ,
                        Z,
                        LDZ,
                        INFO);
}

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
                    const INTEGER    *LDY)
{
    DEBUG_LAPACK_STUB("zlahr2");
    LAPACK_IMPL(zlahr2)(N,
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
                    const INTEGER    *LDY)
{
    DEBUG_LAPACK_STUB("zlahrd");
    LAPACK_IMPL(zlahrd)(N,
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
                    DOUBLE_COMPLEX           *C)
{
    DEBUG_LAPACK_STUB("zlaic1");
    LAPACK_IMPL(zlaic1)(JOB,
                        J,
                        X,
                        SEST,
                        W,
                        GAMMA,
                        SESTPR,
                        S,
                        C);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zlals0");
    LAPACK_IMPL(zlals0)(ICOMPQ,
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
                        RWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zlalsa");
    LAPACK_IMPL(zlalsa)(ICOMPQ,
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
                        RWORK,
                        IWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zlalsd");
    LAPACK_IMPL(zlalsd)(UPLO,
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
                        RWORK,
                        IWORK,
                        INFO);
}

//-- zlangb --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlangb)(const char               *NORM,
                    const INTEGER            *N,
                    const INTEGER            *KL,
                    const INTEGER            *KU,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE                   *WORK)
{
    DEBUG_LAPACK_STUB("zlangb");
    return LAPACK_IMPL(zlangb)(NORM,
                               N,
                               KL,
                               KU,
                               AB,
                               LDAB,
                               WORK);
}

//-- zlange --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlange)(const char               *NORM,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *WORK)
{
    DEBUG_LAPACK_STUB("zlange");
    return LAPACK_IMPL(zlange)(NORM,
                               M,
                               N,
                               A,
                               LDA,
                               WORK);
}

//-- zlangt --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlangt)(const char               *NORM,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *DL,
                    const DOUBLE_COMPLEX     *D,
                    const DOUBLE_COMPLEX     *DU)
{
    DEBUG_LAPACK_STUB("zlangt");
    return LAPACK_IMPL(zlangt)(NORM,
                               N,
                               DL,
                               D,
                               DU);
}

//-- zlanhb --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlanhb)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE                   *WORK)
{
    DEBUG_LAPACK_STUB("zlanhb");
    return LAPACK_IMPL(zlanhb)(NORM,
                               UPLO,
                               N,
                               K,
                               AB,
                               LDAB,
                               WORK);
}

//-- zlanhe --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlanhe)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *WORK)
{
    DEBUG_LAPACK_STUB("zlanhe");
    return LAPACK_IMPL(zlanhe)(NORM,
                               UPLO,
                               N,
                               A,
                               LDA,
                               WORK);
}

//-- zlanhf --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlanhf)(const char               *NORM,
                    const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    DOUBLE                   *WORK)
{
    DEBUG_LAPACK_STUB("zlanhf");
    return LAPACK_IMPL(zlanhf)(NORM,
                               TRANSR,
                               UPLO,
                               N,
                               A,
                               WORK);
}

//-- zlanhp --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlanhp)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE                   *WORK)
{
    DEBUG_LAPACK_STUB("zlanhp");
    return LAPACK_IMPL(zlanhp)(NORM,
                               UPLO,
                               N,
                               AP,
                               WORK);
}

//-- zlanhs --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlanhs)(const char               *NORM,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *WORK)
{
    DEBUG_LAPACK_STUB("zlanhs");
    return LAPACK_IMPL(zlanhs)(NORM,
                               N,
                               A,
                               LDA,
                               WORK);
}

//-- zlanht --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlanht)(const char               *NORM,
                    const INTEGER            *N,
                    const DOUBLE             *D,
                    const DOUBLE_COMPLEX     *E)
{
    DEBUG_LAPACK_STUB("zlanht");
    return LAPACK_IMPL(zlanht)(NORM,
                               N,
                               D,
                               E);
}

//-- zlansb --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlansb)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE                   *WORK)
{
    DEBUG_LAPACK_STUB("zlansb");
    return LAPACK_IMPL(zlansb)(NORM,
                               UPLO,
                               N,
                               K,
                               AB,
                               LDAB,
                               WORK);
}

//-- zlansp --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlansp)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE                   *WORK)
{
    DEBUG_LAPACK_STUB("zlansp");
    return LAPACK_IMPL(zlansp)(NORM,
                               UPLO,
                               N,
                               AP,
                               WORK);
}

//-- zlansy --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlansy)(const char               *NORM,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *WORK)
{
    DEBUG_LAPACK_STUB("zlansy");
    return LAPACK_IMPL(zlansy)(NORM,
                               UPLO,
                               N,
                               A,
                               LDA,
                               WORK);
}

//-- zlantb --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlantb)(const char               *NORM,
                    const char               *UPLO,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    const DOUBLE_COMPLEX     *AB,
                    const INTEGER            *LDAB,
                    DOUBLE                   *WORK)
{
    DEBUG_LAPACK_STUB("zlantb");
    return LAPACK_IMPL(zlantb)(NORM,
                               UPLO,
                               DIAG,
                               N,
                               K,
                               AB,
                               LDAB,
                               WORK);
}

//-- zlantp --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlantp)(const char               *NORM,
                    const char               *UPLO,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE                   *WORK)
{
    DEBUG_LAPACK_STUB("zlantp");
    return LAPACK_IMPL(zlantp)(NORM,
                               UPLO,
                               DIAG,
                               N,
                               AP,
                               WORK);
}

//-- zlantr --------------------------------------------------------------------
DOUBLE
LAPACK_DECL(zlantr)(const char               *NORM,
                    const char               *UPLO,
                    const char               *DIAG,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *WORK)
{
    DEBUG_LAPACK_STUB("zlantr");
    return LAPACK_IMPL(zlantr)(NORM,
                               UPLO,
                               DIAG,
                               M,
                               N,
                               A,
                               LDA,
                               WORK);
}

//-- zlapll --------------------------------------------------------------------
void
LAPACK_DECL(zlapll)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *INCX,
                    DOUBLE_COMPLEX   *Y,
                    const INTEGER    *INCY,
                    DOUBLE           *SSMIN)
{
    DEBUG_LAPACK_STUB("zlapll");
    LAPACK_IMPL(zlapll)(N,
                        X,
                        INCX,
                        Y,
                        INCY,
                        SSMIN);
}

//-- zlapmr --------------------------------------------------------------------
void
LAPACK_DECL(zlapmr)(const LOGICAL    *FORWRD,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *LDX,
                    INTEGER          *K)
{
    DEBUG_LAPACK_STUB("zlapmr");
    LAPACK_IMPL(zlapmr)(FORWRD,
                        M,
                        N,
                        X,
                        LDX,
                        K);
}

//-- zlapmt --------------------------------------------------------------------
void
LAPACK_DECL(zlapmt)(const LOGICAL    *FORWRD,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *LDX,
                    INTEGER          *K)
{
    DEBUG_LAPACK_STUB("zlapmt");
    LAPACK_IMPL(zlapmt)(FORWRD,
                        M,
                        N,
                        X,
                        LDX,
                        K);
}

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
                    char             *EQUED)
{
    DEBUG_LAPACK_STUB("zlaqgb");
    LAPACK_IMPL(zlaqgb)(M,
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
                    char             *EQUED)
{
    DEBUG_LAPACK_STUB("zlaqge");
    LAPACK_IMPL(zlaqge)(M,
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
                    char             *EQUED)
{
    DEBUG_LAPACK_STUB("zlaqhb");
    LAPACK_IMPL(zlaqhb)(UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        S,
                        SCOND,
                        AMAX,
                        EQUED);
}

//-- zlaqhe --------------------------------------------------------------------
void
LAPACK_DECL(zlaqhe)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED)
{
    DEBUG_LAPACK_STUB("zlaqhe");
    LAPACK_IMPL(zlaqhe)(UPLO,
                        N,
                        A,
                        LDA,
                        S,
                        SCOND,
                        AMAX,
                        EQUED);
}

//-- zlaqhp --------------------------------------------------------------------
void
LAPACK_DECL(zlaqhp)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED)
{
    DEBUG_LAPACK_STUB("zlaqhp");
    LAPACK_IMPL(zlaqhp)(UPLO,
                        N,
                        AP,
                        S,
                        SCOND,
                        AMAX,
                        EQUED);
}

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
                    DOUBLE_COMPLEX   *WORK)
{
    DEBUG_LAPACK_STUB("zlaqp2");
    LAPACK_IMPL(zlaqp2)(M,
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
                    const INTEGER    *LDF)
{
    DEBUG_LAPACK_STUB("zlaqps");
    LAPACK_IMPL(zlaqps)(M,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zlaqr0");
    LAPACK_IMPL(zlaqr0)(WANTT,
                        WANTZ,
                        N,
                        ILO,
                        IHI,
                        H,
                        LDH,
                        W,
                        ILOZ,
                        IHIZ,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        INFO);
}

//-- zlaqr1 --------------------------------------------------------------------
void
LAPACK_DECL(zlaqr1)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *H,
                    const INTEGER            *LDH,
                    const DOUBLE_COMPLEX     *S1,
                    const DOUBLE_COMPLEX     *S2,
                    DOUBLE_COMPLEX           *V)
{
    DEBUG_LAPACK_STUB("zlaqr1");
    LAPACK_IMPL(zlaqr1)(N,
                        H,
                        LDH,
                        S1,
                        S2,
                        V);
}

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
                    const INTEGER    *LWORK)
{
    DEBUG_LAPACK_STUB("zlaqr2");
    LAPACK_IMPL(zlaqr2)(WANTT,
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
                        SH,
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
                    const INTEGER    *LWORK)
{
    DEBUG_LAPACK_STUB("zlaqr3");
    LAPACK_IMPL(zlaqr3)(WANTT,
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
                        SH,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zlaqr4");
    LAPACK_IMPL(zlaqr4)(WANTT,
                        WANTZ,
                        N,
                        ILO,
                        IHI,
                        H,
                        LDH,
                        W,
                        ILOZ,
                        IHIZ,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    const INTEGER    *LDWH)
{
    DEBUG_LAPACK_STUB("zlaqr5");
    LAPACK_IMPL(zlaqr5)(WANTT,
                        WANTZ,
                        KACC22,
                        N,
                        KTOP,
                        KBOT,
                        NSHFTS,
                        S,
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
                    char             *EQUED)
{
    DEBUG_LAPACK_STUB("zlaqsb");
    LAPACK_IMPL(zlaqsb)(UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        S,
                        SCOND,
                        AMAX,
                        EQUED);
}

//-- zlaqsp --------------------------------------------------------------------
void
LAPACK_DECL(zlaqsp)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED)
{
    DEBUG_LAPACK_STUB("zlaqsp");
    LAPACK_IMPL(zlaqsp)(UPLO,
                        N,
                        AP,
                        S,
                        SCOND,
                        AMAX,
                        EQUED);
}

//-- zlaqsy --------------------------------------------------------------------
void
LAPACK_DECL(zlaqsy)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *S,
                    const DOUBLE     *SCOND,
                    const DOUBLE     *AMAX,
                    char             *EQUED)
{
    DEBUG_LAPACK_STUB("zlaqsy");
    LAPACK_IMPL(zlaqsy)(UPLO,
                        N,
                        A,
                        LDA,
                        S,
                        SCOND,
                        AMAX,
                        EQUED);
}

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
                    DOUBLE           *WORK)
{
    DEBUG_LAPACK_STUB("zlar1v");
    LAPACK_IMPL(zlar1v)(N,
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

//-- zlar2v --------------------------------------------------------------------
void
LAPACK_DECL(zlar2v)(const INTEGER            *N,
                    DOUBLE_COMPLEX           *X,
                    DOUBLE_COMPLEX           *Y,
                    DOUBLE_COMPLEX           *Z,
                    const INTEGER            *INCX,
                    const DOUBLE             *C,
                    const DOUBLE_COMPLEX     *S,
                    const INTEGER            *INCC)
{
    DEBUG_LAPACK_STUB("zlar2v");
    LAPACK_IMPL(zlar2v)(N,
                        X,
                        Y,
                        Z,
                        INCX,
                        C,
                        S,
                        INCC);
}

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
                    DOUBLE                   *RWORK)
{
    DEBUG_LAPACK_STUB("zlarcm");
    LAPACK_IMPL(zlarcm)(M,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        C,
                        LDC,
                        RWORK);
}

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
                   DOUBLE_COMPLEX           *WORK)
{
    DEBUG_LAPACK_STUB("zlarf");
    LAPACK_IMPL(zlarf)(SIDE,
                       M,
                       N,
                       V,
                       INCV,
                       TAU,
                       C,
                       LDC,
                       WORK);
}

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
                    const INTEGER            *LDWORK)
{
    DEBUG_LAPACK_STUB("zlarfb");
    LAPACK_IMPL(zlarfb)(SIDE,
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

//-- zlarfg --------------------------------------------------------------------
void
LAPACK_DECL(zlarfg)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *ALPHA,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *INCX,
                    DOUBLE_COMPLEX   *TAU)
{
    DEBUG_LAPACK_STUB("zlarfg");
    LAPACK_IMPL(zlarfg)(N,
                        ALPHA,
                        X,
                        INCX,
                        TAU);
}

//-- zlarfgp -------------------------------------------------------------------
void
LAPACK_DECL(zlarfgp)(const INTEGER    *N,
                     DOUBLE_COMPLEX   *ALPHA,
                     DOUBLE_COMPLEX   *X,
                     const INTEGER    *INCX,
                     DOUBLE_COMPLEX   *TAU)
{
    DEBUG_LAPACK_STUB("zlarfgp");
    LAPACK_IMPL(zlarfgp)(N,
                         ALPHA,
                         X,
                         INCX,
                         TAU);
}

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
                    const INTEGER            *LDT)
{
    DEBUG_LAPACK_STUB("zlarft");
    LAPACK_IMPL(zlarft)(DIRECT,
                        STOREV,
                        N,
                        K,
                        V,
                        LDV,
                        TAU,
                        T,
                        LDT);
}

//-- zlarfx --------------------------------------------------------------------
void
LAPACK_DECL(zlarfx)(const char               *SIDE,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *V,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK)
{
    DEBUG_LAPACK_STUB("zlarfx");
    LAPACK_IMPL(zlarfx)(SIDE,
                        M,
                        N,
                        V,
                        TAU,
                        C,
                        LDC,
                        WORK);
}

//-- zlargv --------------------------------------------------------------------
void
LAPACK_DECL(zlargv)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *X,
                    const INTEGER    *INCX,
                    DOUBLE_COMPLEX   *Y,
                    const INTEGER    *INCY,
                    DOUBLE           *C,
                    const INTEGER    *INCC)
{
    DEBUG_LAPACK_STUB("zlargv");
    LAPACK_IMPL(zlargv)(N,
                        X,
                        INCX,
                        Y,
                        INCY,
                        C,
                        INCC);
}

//-- zlarnv --------------------------------------------------------------------
void
LAPACK_DECL(zlarnv)(const INTEGER    *IDIST,
                    INTEGER          *ISEED,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *X)
{
    DEBUG_LAPACK_STUB("zlarnv");
    LAPACK_IMPL(zlarnv)(IDIST,
                        ISEED,
                        N,
                        X);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zlarrv");
    LAPACK_IMPL(zlarrv)(N,
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

//-- zlarscl2 ------------------------------------------------------------------
void
LAPACK_DECL(zlarscl2)(const INTEGER        *M,
                      const INTEGER        *N,
                      const DOUBLE         *D,
                      DOUBLE_COMPLEX       *X,
                      const INTEGER        *LDX)
{
    DEBUG_LAPACK_STUB("zlarscl2");
    LAPACK_IMPL(zlarscl2)(M,
                          N,
                          D,
                          X,
                          LDX);
}

//-- zlartg --------------------------------------------------------------------
void
LAPACK_DECL(zlartg)(const DOUBLE_COMPLEX     *F,
                    const DOUBLE_COMPLEX     *G,
                    DOUBLE                   *CS,
                    DOUBLE_COMPLEX           *SN,
                    DOUBLE_COMPLEX           *R)
{
    DEBUG_LAPACK_STUB("zlartg");
    LAPACK_IMPL(zlartg)(F,
                        G,
                        CS,
                        SN,
                        R);
}

//-- zlartv --------------------------------------------------------------------
void
LAPACK_DECL(zlartv)(const INTEGER            *N,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *INCX,
                    DOUBLE_COMPLEX           *Y,
                    const INTEGER            *INCY,
                    const DOUBLE             *C,
                    const DOUBLE_COMPLEX     *S,
                    const INTEGER            *INCC)
{
    DEBUG_LAPACK_STUB("zlartv");
    LAPACK_IMPL(zlartv)(N,
                        X,
                        INCX,
                        Y,
                        INCY,
                        C,
                        S,
                        INCC);
}

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
                   DOUBLE_COMPLEX           *WORK)
{
    DEBUG_LAPACK_STUB("zlarz");
    LAPACK_IMPL(zlarz)(SIDE,
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
                    const INTEGER            *LDWORK)
{
    DEBUG_LAPACK_STUB("zlarzb");
    LAPACK_IMPL(zlarzb)(SIDE,
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
                    const INTEGER            *LDT)
{
    DEBUG_LAPACK_STUB("zlarzt");
    LAPACK_IMPL(zlarzt)(DIRECT,
                        STOREV,
                        N,
                        K,
                        V,
                        LDV,
                        TAU,
                        T,
                        LDT);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zlascl");
    LAPACK_IMPL(zlascl)(TYPE,
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

//-- zlascl2 -------------------------------------------------------------------
void
LAPACK_DECL(zlascl2)(const INTEGER    *M,
                     const INTEGER    *N,
                     const DOUBLE     *D,
                     DOUBLE_COMPLEX   *X,
                     const INTEGER    *LDX)
{
    DEBUG_LAPACK_STUB("zlascl2");
    LAPACK_IMPL(zlascl2)(M,
                         N,
                         D,
                         X,
                         LDX);
}

//-- zlaset --------------------------------------------------------------------
void
LAPACK_DECL(zlaset)(const char               *UPLO,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *ALPHA,
                    const DOUBLE_COMPLEX     *BETA,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA)
{
    DEBUG_LAPACK_STUB("zlaset");
    LAPACK_IMPL(zlaset)(UPLO,
                        M,
                        N,
                        ALPHA,
                        BETA,
                        A,
                        LDA);
}

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
                   const INTEGER        *LDA)
{
    DEBUG_LAPACK_STUB("zlasr");
    LAPACK_IMPL(zlasr)(SIDE,
                       PIVOT,
                       DIRECT,
                       M,
                       N,
                       C,
                       S,
                       A,
                       LDA);
}

//-- zlassq --------------------------------------------------------------------
void
LAPACK_DECL(zlassq)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *X,
                    const INTEGER            *INCX,
                    DOUBLE                   *SCALE,
                    DOUBLE                   *SUMSQ)
{
    DEBUG_LAPACK_STUB("zlassq");
    LAPACK_IMPL(zlassq)(N,
                        X,
                        INCX,
                        SCALE,
                        SUMSQ);
}

//-- zlaswp --------------------------------------------------------------------
void
LAPACK_DECL(zlaswp)(const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const INTEGER    *K1,
                    const INTEGER    *K2,
                    const INTEGER    *IPIV,
                    const INTEGER    *INCX)
{
    DEBUG_LAPACK_STUB("zlaswp");
    LAPACK_IMPL(zlaswp)(N,
                        A,
                        LDA,
                        K1,
                        K2,
                        IPIV,
                        INCX);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zlasyf");
    LAPACK_IMPL(zlasyf)(UPLO,
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

//-- zlat2c --------------------------------------------------------------------
void
LAPACK_DECL(zlat2c)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    FLOAT_COMPLEX            *SA,
                    const INTEGER            *LDSA,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zlat2c");
    LAPACK_IMPL(zlat2c)(UPLO,
                        N,
                        A,
                        LDA,
                        SA,
                        LDSA,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zlatbs");
    LAPACK_IMPL(zlatbs)(UPLO,
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
                    const INTEGER            *JPIV)
{
    DEBUG_LAPACK_STUB("zlatdf");
    LAPACK_IMPL(zlatdf)(IJOB,
                        N,
                        Z,
                        LDZ,
                        RHS,
                        RDSUM,
                        RDSCAL,
                        IPIV,
                        JPIV);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zlatps");
    LAPACK_IMPL(zlatps)(UPLO,
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
                    const INTEGER    *LDW)
{
    DEBUG_LAPACK_STUB("zlatrd");
    LAPACK_IMPL(zlatrd)(UPLO,
                        N,
                        NB,
                        A,
                        LDA,
                        E,
                        TAU,
                        W,
                        LDW);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zlatrs");
    LAPACK_IMPL(zlatrs)(UPLO,
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

//-- zlatrz --------------------------------------------------------------------
void
LAPACK_DECL(zlatrz)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *L,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK)
{
    DEBUG_LAPACK_STUB("zlatrz");
    LAPACK_IMPL(zlatrz)(M,
                        N,
                        L,
                        A,
                        LDA,
                        TAU,
                        WORK);
}

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
                    DOUBLE_COMPLEX           *WORK)
{
    DEBUG_LAPACK_STUB("zlatzm");
    LAPACK_IMPL(zlatzm)(SIDE,
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

//-- zlauu2 --------------------------------------------------------------------
void
LAPACK_DECL(zlauu2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zlauu2");
    LAPACK_IMPL(zlauu2)(UPLO,
                        N,
                        A,
                        LDA,
                        INFO);
}

//-- zlauum --------------------------------------------------------------------
void
LAPACK_DECL(zlauum)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zlauum");
    LAPACK_IMPL(zlauum)(UPLO,
                        N,
                        A,
                        LDA,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zpbcon");
    LAPACK_IMPL(zpbcon)(UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        ANORM,
                        RCOND,
                        WORK,
                        RWORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zpbequ");
    LAPACK_IMPL(zpbequ)(UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        S,
                        SCOND,
                        AMAX,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zpbrfs");
    LAPACK_IMPL(zpbrfs)(UPLO,
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
                        RWORK,
                        INFO);
}

//-- zpbstf --------------------------------------------------------------------
void
LAPACK_DECL(zpbstf)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zpbstf");
    LAPACK_IMPL(zpbstf)(UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        INFO);
}

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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zpbsv");
    LAPACK_IMPL(zpbsv)(UPLO,
                       N,
                       KD,
                       NRHS,
                       AB,
                       LDAB,
                       B,
                       LDB,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zpbsvx");
    LAPACK_IMPL(zpbsvx)(FACT,
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
                        RWORK,
                        INFO);
}

//-- zpbtf2 --------------------------------------------------------------------
void
LAPACK_DECL(zpbtf2)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zpbtf2");
    LAPACK_IMPL(zpbtf2)(UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        INFO);
}

//-- zpbtrf --------------------------------------------------------------------
void
LAPACK_DECL(zpbtrf)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *KD,
                    DOUBLE_COMPLEX   *AB,
                    const INTEGER    *LDAB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zpbtrf");
    LAPACK_IMPL(zpbtrf)(UPLO,
                        N,
                        KD,
                        AB,
                        LDAB,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zpbtrs");
    LAPACK_IMPL(zpbtrs)(UPLO,
                        N,
                        KD,
                        NRHS,
                        AB,
                        LDAB,
                        B,
                        LDB,
                        INFO);
}

//-- zpftrf --------------------------------------------------------------------
void
LAPACK_DECL(zpftrf)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zpftrf");
    LAPACK_IMPL(zpftrf)(TRANSR,
                        UPLO,
                        N,
                        A,
                        INFO);
}

//-- zpftri --------------------------------------------------------------------
void
LAPACK_DECL(zpftri)(const char       *TRANSR,
                    const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zpftri");
    LAPACK_IMPL(zpftri)(TRANSR,
                        UPLO,
                        N,
                        A,
                        INFO);
}

//-- zpftrs --------------------------------------------------------------------
void
LAPACK_DECL(zpftrs)(const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zpftrs");
    LAPACK_IMPL(zpftrs)(TRANSR,
                        UPLO,
                        N,
                        NRHS,
                        A,
                        B,
                        LDB,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zpocon");
    LAPACK_IMPL(zpocon)(UPLO,
                        N,
                        A,
                        LDA,
                        ANORM,
                        RCOND,
                        WORK,
                        RWORK,
                        INFO);
}

//-- zpoequ --------------------------------------------------------------------
void
LAPACK_DECL(zpoequ)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *S,
                    DOUBLE                   *SCOND,
                    DOUBLE                   *AMAX,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zpoequ");
    LAPACK_IMPL(zpoequ)(N,
                        A,
                        LDA,
                        S,
                        SCOND,
                        AMAX,
                        INFO);
}

//-- zpoequb -------------------------------------------------------------------
void
LAPACK_DECL(zpoequb)(const INTEGER            *N,
                     const DOUBLE_COMPLEX     *A,
                     const INTEGER            *LDA,
                     DOUBLE                   *S,
                     DOUBLE                   *SCOND,
                     DOUBLE                   *AMAX,
                     INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zpoequb");
    LAPACK_IMPL(zpoequb)(N,
                         A,
                         LDA,
                         S,
                         SCOND,
                         AMAX,
                         INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zporfs");
    LAPACK_IMPL(zporfs)(UPLO,
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
                        RWORK,
                        INFO);
}

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
                     INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zporfsx");
    LAPACK_IMPL(zporfsx)(UPLO,
                         EQUED,
                         N,
                         NRHS,
                         A,
                         LDA,
                         AF,
                         LDAF,
                         S,
                         B,
                         LDB,
                         X,
                         LDX,
                         RCOND,
                         BERR,
                         N_ERR_BNDS,
                         ERR_BNDS_NORM,
                         ERR_BNDS_COMP,
                         NPARAMS,
                         PARAMS,
                         WORK,
                         RWORK,
                         INFO);
}
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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zposv");
    LAPACK_IMPL(zposv)(UPLO,
                       N,
                       NRHS,
                       A,
                       LDA,
                       B,
                       LDB,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zposvx");
    LAPACK_IMPL(zposvx)(FACT,
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
                        RWORK,
                        INFO);
}

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
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zposvxx");
    LAPACK_IMPL(zposvxx)(FACT,
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
                         RPVGRW,
                         BERR,
                         N_ERR_BNDS,
                         ERR_BNDS_NORM,
                         ERR_BNDS_COMP,
                         NPARAMS,
                         PARAMS,
                         WORK,
                         RWORK,
                         INFO);
}
*/

//-- zpotf2 --------------------------------------------------------------------
void
LAPACK_DECL(zpotf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zpotf2");
    LAPACK_IMPL(zpotf2)(UPLO,
                        N,
                        A,
                        LDA,
                        INFO);
}

//-- zpotrf --------------------------------------------------------------------
void
LAPACK_DECL(zpotrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zpotrf");
    LAPACK_IMPL(zpotrf)(UPLO,
                        N,
                        A,
                        LDA,
                        INFO);
}

//-- zpotri --------------------------------------------------------------------
void
LAPACK_DECL(zpotri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zpotri");
    LAPACK_IMPL(zpotri)(UPLO,
                        N,
                        A,
                        LDA,
                        INFO);
}

//-- zpotrs --------------------------------------------------------------------
void
LAPACK_DECL(zpotrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zpotrs");
    LAPACK_IMPL(zpotrs)(UPLO,
                        N,
                        NRHS,
                        A,
                        LDA,
                        B,
                        LDB,
                        INFO);
}

//-- zppcon --------------------------------------------------------------------
void
LAPACK_DECL(zppcon)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zppcon");
    LAPACK_IMPL(zppcon)(UPLO,
                        N,
                        AP,
                        ANORM,
                        RCOND,
                        WORK,
                        RWORK,
                        INFO);
}

//-- zppequ --------------------------------------------------------------------
void
LAPACK_DECL(zppequ)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE                   *S,
                    DOUBLE                   *SCOND,
                    DOUBLE                   *AMAX,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zppequ");
    LAPACK_IMPL(zppequ)(UPLO,
                        N,
                        AP,
                        S,
                        SCOND,
                        AMAX,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zpprfs");
    LAPACK_IMPL(zpprfs)(UPLO,
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
                        RWORK,
                        INFO);
}

//-- zppsv ---------------------------------------------------------------------
void
LAPACK_DECL(zppsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *AP,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zppsv");
    LAPACK_IMPL(zppsv)(UPLO,
                       N,
                       NRHS,
                       AP,
                       B,
                       LDB,
                       INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zppsvx");
    LAPACK_IMPL(zppsvx)(FACT,
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
                        RWORK,
                        INFO);
}

//-- zpptrf --------------------------------------------------------------------
void
LAPACK_DECL(zpptrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zpptrf");
    LAPACK_IMPL(zpptrf)(UPLO,
                        N,
                        AP,
                        INFO);
}

//-- zpptri --------------------------------------------------------------------
void
LAPACK_DECL(zpptri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zpptri");
    LAPACK_IMPL(zpptri)(UPLO,
                        N,
                        AP,
                        INFO);
}

//-- zpptrs --------------------------------------------------------------------
void
LAPACK_DECL(zpptrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zpptrs");
    LAPACK_IMPL(zpptrs)(UPLO,
                        N,
                        NRHS,
                        AP,
                        B,
                        LDB,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zpstf2");
    LAPACK_IMPL(zpstf2)(UPLO,
                        N,
                        A,
                        LDA,
                        PIV,
                        RANK,
                        TOL,
                        WORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zpstrf");
    LAPACK_IMPL(zpstrf)(UPLO,
                        N,
                        A,
                        LDA,
                        PIV,
                        RANK,
                        TOL,
                        WORK,
                        INFO);
}

//-- zptcon --------------------------------------------------------------------
void
LAPACK_DECL(zptcon)(const INTEGER            *N,
                    const DOUBLE             *D,
                    const DOUBLE_COMPLEX     *E,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zptcon");
    LAPACK_IMPL(zptcon)(N,
                        D,
                        E,
                        ANORM,
                        RCOND,
                        RWORK,
                        INFO);
}

//-- zpteqr --------------------------------------------------------------------
void
LAPACK_DECL(zpteqr)(const char       *COMPZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zpteqr");
    LAPACK_IMPL(zpteqr)(COMPZ,
                        N,
                        D,
                        E,
                        Z,
                        LDZ,
                        WORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zptrfs");
    LAPACK_IMPL(zptrfs)(UPLO,
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
                        FERR,
                        BERR,
                        WORK,
                        RWORK,
                        INFO);
}

//-- zptsv ---------------------------------------------------------------------
void
LAPACK_DECL(zptsv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *D,
                   DOUBLE_COMPLEX       *E,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zptsv");
    LAPACK_IMPL(zptsv)(N,
                       NRHS,
                       D,
                       E,
                       B,
                       LDB,
                       INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zptsvx");
    LAPACK_IMPL(zptsvx)(FACT,
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
                        RWORK,
                        INFO);
}

//-- zpttrf --------------------------------------------------------------------
void
LAPACK_DECL(zpttrf)(const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE_COMPLEX   *E,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zpttrf");
    LAPACK_IMPL(zpttrf)(N,
                        D,
                        E,
                        INFO);
}

//-- zpttrs --------------------------------------------------------------------
void
LAPACK_DECL(zpttrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE             *D,
                    const DOUBLE_COMPLEX     *E,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zpttrs");
    LAPACK_IMPL(zpttrs)(UPLO,
                        N,
                        NRHS,
                        D,
                        E,
                        B,
                        LDB,
                        INFO);
}

//-- zptts2 --------------------------------------------------------------------
void
LAPACK_DECL(zptts2)(const INTEGER            *IUPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE             *D,
                    const DOUBLE_COMPLEX     *E,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB)
{
    DEBUG_LAPACK_STUB("zptts2");
    LAPACK_IMPL(zptts2)(IUPLO,
                        N,
                        NRHS,
                        D,
                        E,
                        B,
                        LDB);
}

//-- zrot ----------------------------------------------------------------------
void
LAPACK_DECL(zrot)(const INTEGER            *N,
                  DOUBLE_COMPLEX           *CX,
                  const INTEGER            *INCX,
                  DOUBLE_COMPLEX           *CY,
                  const INTEGER            *INCY,
                  const DOUBLE             *C,
                  const DOUBLE_COMPLEX     *S)
{
    DEBUG_LAPACK_STUB("zrot");
    LAPACK_IMPL(zrot)(N,
                      CX,
                      INCX,
                      CY,
                      INCY,
                      C,
                      S);
}

//-- zspcon --------------------------------------------------------------------
void
LAPACK_DECL(zspcon)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    const INTEGER            *IPIV,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zspcon");
    LAPACK_IMPL(zspcon)(UPLO,
                        N,
                        AP,
                        IPIV,
                        ANORM,
                        RCOND,
                        WORK,
                        INFO);
}

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
                   const INTEGER            *INCY)
{
    DEBUG_LAPACK_STUB("zspmv");
    LAPACK_IMPL(zspmv)(UPLO,
                       N,
                       ALPHA,
                       AP,
                       X,
                       INCX,
                       BETA,
                       Y,
                       INCY);
}

//-- zspr ----------------------------------------------------------------------
void
LAPACK_DECL(zspr)(const char               *UPLO,
                  const INTEGER            *N,
                  const DOUBLE_COMPLEX     *ALPHA,
                  const DOUBLE_COMPLEX     *X,
                  const INTEGER            *INCX,
                  DOUBLE_COMPLEX           *AP)
{
    DEBUG_LAPACK_STUB("zspr");
    LAPACK_IMPL(zspr)(UPLO,
                      N,
                      ALPHA,
                      X,
                      INCX,
                      AP);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zsprfs");
    LAPACK_IMPL(zsprfs)(UPLO,
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
                        RWORK,
                        INFO);
}

//-- zspsv ---------------------------------------------------------------------
void
LAPACK_DECL(zspsv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *AP,
                   INTEGER              *IPIV,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zspsv");
    LAPACK_IMPL(zspsv)(UPLO,
                       N,
                       NRHS,
                       AP,
                       IPIV,
                       B,
                       LDB,
                       INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zspsvx");
    LAPACK_IMPL(zspsvx)(FACT,
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
                        RWORK,
                        INFO);
}

//-- zsptrf --------------------------------------------------------------------
void
LAPACK_DECL(zsptrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zsptrf");
    LAPACK_IMPL(zsptrf)(UPLO,
                        N,
                        AP,
                        IPIV,
                        INFO);
}

//-- zsptri --------------------------------------------------------------------
void
LAPACK_DECL(zsptri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    const INTEGER    *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zsptri");
    LAPACK_IMPL(zsptri)(UPLO,
                        N,
                        AP,
                        IPIV,
                        WORK,
                        INFO);
}

//-- zsptrs --------------------------------------------------------------------
void
LAPACK_DECL(zsptrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *AP,
                    const INTEGER            *IPIV,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zsptrs");
    LAPACK_IMPL(zsptrs)(UPLO,
                        N,
                        NRHS,
                        AP,
                        IPIV,
                        B,
                        LDB,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zstedc");
    LAPACK_IMPL(zstedc)(COMPZ,
                        N,
                        D,
                        E,
                        Z,
                        LDZ,
                        WORK,
                        LWORK,
                        RWORK,
                        LRWORK,
                        IWORK,
                        LIWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zstegr");
    LAPACK_IMPL(zstegr)(JOBZ,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zstein");
    LAPACK_IMPL(zstein)(N,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zstemr");
    LAPACK_IMPL(zstemr)(JOBZ,
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

//-- zsteqr --------------------------------------------------------------------
void
LAPACK_DECL(zsteqr)(const char       *COMPZ,
                    const INTEGER    *N,
                    DOUBLE           *D,
                    DOUBLE           *E,
                    DOUBLE_COMPLEX   *Z,
                    const INTEGER    *LDZ,
                    DOUBLE           *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zsteqr");
    LAPACK_IMPL(zsteqr)(COMPZ,
                        N,
                        D,
                        E,
                        Z,
                        LDZ,
                        WORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zsycon");
    LAPACK_IMPL(zsycon)(UPLO,
                        N,
                        A,
                        LDA,
                        IPIV,
                        ANORM,
                        RCOND,
                        WORK,
                        INFO);
}

//-- zsyconv -------------------------------------------------------------------
void
LAPACK_DECL(zsyconv)(const char       *UPLO,
                     const char       *WAY,
                     const INTEGER    *N,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE_COMPLEX   *WORK,
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zsyconv");
    LAPACK_IMPL(zsyconv)(UPLO,
                         WAY,
                         N,
                         A,
                         LDA,
                         IPIV,
                         WORK,
                         INFO);
}

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
                     INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zsyequb");
    LAPACK_IMPL(zsyequb)(UPLO,
                         N,
                         A,
                         LDA,
                         S,
                         SCOND,
                         AMAX,
                         WORK,
                         INFO);
}

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
                   const INTEGER            *INCY)
{
    DEBUG_LAPACK_STUB("zsymv");
    LAPACK_IMPL(zsymv)(UPLO,
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

//-- zsyr ----------------------------------------------------------------------
void
LAPACK_DECL(zsyr)(const char               *UPLO,
                  const INTEGER            *N,
                  const DOUBLE_COMPLEX     *ALPHA,
                  const DOUBLE_COMPLEX     *X,
                  const INTEGER            *INCX,
                  DOUBLE_COMPLEX           *A,
                  const INTEGER            *LDA)
{
    DEBUG_LAPACK_STUB("zsyr");
    LAPACK_IMPL(zsyr)(UPLO,
                      N,
                      ALPHA,
                      X,
                      INCX,
                      A,
                      LDA);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zsyrfs");
    LAPACK_IMPL(zsyrfs)(UPLO,
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
                        RWORK,
                        INFO);
}

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
                     INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zsyrfsx");
    LAPACK_IMPL(zsyrfsx)(UPLO,
                         EQUED,
                         N,
                         NRHS,
                         A,
                         LDA,
                         AF,
                         LDAF,
                         IPIV,
                         S,
                         B,
                         LDB,
                         X,
                         LDX,
                         RCOND,
                         BERR,
                         N_ERR_BNDS,
                         ERR_BNDS_NORM,
                         ERR_BNDS_COMP,
                         NPARAMS,
                         PARAMS,
                         WORK,
                         RWORK,
                         INFO);
}
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
                   INTEGER              *INFO)
{
    DEBUG_LAPACK_STUB("zsysv");
    LAPACK_IMPL(zsysv)(UPLO,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zsysvx");
    LAPACK_IMPL(zsysvx)(FACT,
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
                        RWORK,
                        INFO);
}

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
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zsysvxx");
    LAPACK_IMPL(zsysvxx)(FACT,
                         UPLO,
                         N,
                         NRHS,
                         A,
                         LDA,
                         AF,
                         LDAF,
                         IPIV,
                         EQUED,
                         S,
                         B,
                         LDB,
                         X,
                         LDX,
                         RCOND,
                         RPVGRW,
                         BERR,
                         N_ERR_BNDS,
                         ERR_BNDS_NORM,
                         ERR_BNDS_COMP,
                         NPARAMS,
                         PARAMS,
                         WORK,
                         RWORK,
                         INFO);
}
*/

//-- zsyswapr ------------------------------------------------------------------
void
LAPACK_DECL(zsyswapr)(const char       *UPLO,
                      const INTEGER    *N,
                      DOUBLE_COMPLEX   *A,
                      const INTEGER    *LDA,
                      const INTEGER    *I1,
                      const INTEGER    *I2)
{
    DEBUG_LAPACK_STUB("zsyswapr");
    LAPACK_IMPL(zsyswapr)(UPLO,
                          N,
                          A,
                          LDA,
                          I1,
                          I2);
}

//-- zsytf2 --------------------------------------------------------------------
void
LAPACK_DECL(zsytf2)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zsytf2");
    LAPACK_IMPL(zsytf2)(UPLO,
                        N,
                        A,
                        LDA,
                        IPIV,
                        INFO);
}

//-- zsytrf --------------------------------------------------------------------
void
LAPACK_DECL(zsytrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zsytrf");
    LAPACK_IMPL(zsytrf)(UPLO,
                        N,
                        A,
                        LDA,
                        IPIV,
                        WORK,
                        LWORK,
                        INFO);
}

//-- zsytri --------------------------------------------------------------------
void
LAPACK_DECL(zsytri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE_COMPLEX   *WORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zsytri");
    LAPACK_IMPL(zsytri)(UPLO,
                        N,
                        A,
                        LDA,
                        IPIV,
                        WORK,
                        INFO);
}

//-- zsytri2 -------------------------------------------------------------------
void
LAPACK_DECL(zsytri2)(const char       *UPLO,
                     const INTEGER    *N,
                     DOUBLE_COMPLEX   *A,
                     const INTEGER    *LDA,
                     const INTEGER    *IPIV,
                     DOUBLE_COMPLEX   *WORK,
                     const INTEGER    *LWORK,
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zsytri2");
    LAPACK_IMPL(zsytri2)(UPLO,
                         N,
                         A,
                         LDA,
                         IPIV,
                         WORK,
                         LWORK,
                         INFO);
}

//-- zsytri2x ------------------------------------------------------------------
void
LAPACK_DECL(zsytri2x)(const char       *UPLO,
                      const INTEGER    *N,
                      DOUBLE_COMPLEX   *A,
                      const INTEGER    *LDA,
                      const INTEGER    *IPIV,
                      DOUBLE_COMPLEX   *WORK,
                      const INTEGER    *NB,
                      INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zsytri2x");
    LAPACK_IMPL(zsytri2x)(UPLO,
                          N,
                          A,
                          LDA,
                          IPIV,
                          WORK,
                          NB,
                          INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zsytrs");
    LAPACK_IMPL(zsytrs)(UPLO,
                        N,
                        NRHS,
                        A,
                        LDA,
                        IPIV,
                        B,
                        LDB,
                        INFO);
}

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
                     INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zsytrs2");
    LAPACK_IMPL(zsytrs2)(UPLO,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztbcon");
    LAPACK_IMPL(ztbcon)(NORM,
                        UPLO,
                        DIAG,
                        N,
                        KD,
                        AB,
                        LDAB,
                        RCOND,
                        WORK,
                        RWORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztbrfs");
    LAPACK_IMPL(ztbrfs)(UPLO,
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
                        RWORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztbtrs");
    LAPACK_IMPL(ztbtrs)(UPLO,
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
                   const INTEGER            *LDB)
{
    DEBUG_LAPACK_STUB("ztfsm");
    LAPACK_IMPL(ztfsm)(TRANSR,
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

//-- ztftri --------------------------------------------------------------------
void
LAPACK_DECL(ztftri)(const char       *TRANSR,
                    const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("ztftri");
    LAPACK_IMPL(ztftri)(TRANSR,
                        UPLO,
                        DIAG,
                        N,
                        A,
                        INFO);
}

//-- ztfttp --------------------------------------------------------------------
void
LAPACK_DECL(ztfttp)(const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *ARF,
                    DOUBLE_COMPLEX           *AP,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztfttp");
    LAPACK_IMPL(ztfttp)(TRANSR,
                        UPLO,
                        N,
                        ARF,
                        AP,
                        INFO);
}

//-- ztfttr --------------------------------------------------------------------
void
LAPACK_DECL(ztfttr)(const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *ARF,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztfttr");
    LAPACK_IMPL(ztfttr)(TRANSR,
                        UPLO,
                        N,
                        ARF,
                        A,
                        LDA,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztgevc");
    LAPACK_IMPL(ztgevc)(SIDE,
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
                        RWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("ztgex2");
    LAPACK_IMPL(ztgex2)(WANTQ,
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
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("ztgexc");
    LAPACK_IMPL(ztgexc)(WANTQ,
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
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("ztgsen");
    LAPACK_IMPL(ztgsen)(IJOB,
                        WANTQ,
                        WANTZ,
                        SELECT,
                        N,
                        A,
                        LDA,
                        B,
                        LDB,
                        ALPHA,
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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("ztgsja");
    LAPACK_IMPL(ztgsja)(JOBU,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztgsna");
    LAPACK_IMPL(ztgsna)(JOB,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztgsy2");
    LAPACK_IMPL(ztgsy2)(TRANS,
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
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztgsyl");
    LAPACK_IMPL(ztgsyl)(TRANS,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztpcon");
    LAPACK_IMPL(ztpcon)(NORM,
                        UPLO,
                        DIAG,
                        N,
                        AP,
                        RCOND,
                        WORK,
                        RWORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztprfs");
    LAPACK_IMPL(ztprfs)(UPLO,
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
                        RWORK,
                        INFO);
}

//-- ztptri --------------------------------------------------------------------
void
LAPACK_DECL(ztptri)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *AP,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("ztptri");
    LAPACK_IMPL(ztptri)(UPLO,
                        DIAG,
                        N,
                        AP,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztptrs");
    LAPACK_IMPL(ztptrs)(UPLO,
                        TRANS,
                        DIAG,
                        N,
                        NRHS,
                        AP,
                        B,
                        LDB,
                        INFO);
}

//-- ztpttf --------------------------------------------------------------------
void
LAPACK_DECL(ztpttf)(const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE_COMPLEX           *ARF,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztpttf");
    LAPACK_IMPL(ztpttf)(TRANSR,
                        UPLO,
                        N,
                        AP,
                        ARF,
                        INFO);
}

//-- ztpttr --------------------------------------------------------------------
void
LAPACK_DECL(ztpttr)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztpttr");
    LAPACK_IMPL(ztpttr)(UPLO,
                        N,
                        AP,
                        A,
                        LDA,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztrcon");
    LAPACK_IMPL(ztrcon)(NORM,
                        UPLO,
                        DIAG,
                        N,
                        A,
                        LDA,
                        RCOND,
                        WORK,
                        RWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("ztrevc");
    LAPACK_IMPL(ztrevc)(SIDE,
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
                        RWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("ztrexc");
    LAPACK_IMPL(ztrexc)(COMPQ,
                        N,
                        T,
                        LDT,
                        Q,
                        LDQ,
                        IFST,
                        ILST,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztrrfs");
    LAPACK_IMPL(ztrrfs)(UPLO,
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
                        RWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("ztrsen");
    LAPACK_IMPL(ztrsen)(JOB,
                        COMPQ,
                        SELECT,
                        N,
                        T,
                        LDT,
                        Q,
                        LDQ,
                        W,
                        M,
                        S,
                        SEP,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztrsna");
    LAPACK_IMPL(ztrsna)(JOB,
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
                        RWORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztrsyl");
    LAPACK_IMPL(ztrsyl)(TRANA,
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

//-- ztrti2 --------------------------------------------------------------------
void
LAPACK_DECL(ztrti2)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("ztrti2");
    LAPACK_IMPL(ztrti2)(UPLO,
                        DIAG,
                        N,
                        A,
                        LDA,
                        INFO);
}

//-- ztrtri --------------------------------------------------------------------
void
LAPACK_DECL(ztrtri)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("ztrtri");
    LAPACK_IMPL(ztrtri)(UPLO,
                        DIAG,
                        N,
                        A,
                        LDA,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztrtrs");
    LAPACK_IMPL(ztrtrs)(UPLO,
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

//-- ztrttf --------------------------------------------------------------------
void
LAPACK_DECL(ztrttf)(const char               *TRANSR,
                    const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *ARF,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztrttf");
    LAPACK_IMPL(ztrttf)(TRANSR,
                        UPLO,
                        N,
                        A,
                        LDA,
                        ARF,
                        INFO);
}

//-- ztrttp --------------------------------------------------------------------
void
LAPACK_DECL(ztrttp)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *AP,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("ztrttp");
    LAPACK_IMPL(ztrttp)(UPLO,
                        N,
                        A,
                        LDA,
                        AP,
                        INFO);
}

//-- ztzrqf --------------------------------------------------------------------
void
LAPACK_DECL(ztzrqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("ztzrqf");
    LAPACK_IMPL(ztzrqf)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        INFO);
}

//-- ztzrzf --------------------------------------------------------------------
void
LAPACK_DECL(ztzrzf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("ztzrzf");
    LAPACK_IMPL(ztzrzf)(M,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("zunbdb");
    LAPACK_IMPL(zunbdb)(TRANS,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zuncsd");
    LAPACK_IMPL(zuncsd)(JOBU1,
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
                        RWORK,
                        LRWORK,
                        IWORK,
                        INFO);
}

//-- zung2l --------------------------------------------------------------------
void
LAPACK_DECL(zung2l)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zung2l");
    LAPACK_IMPL(zung2l)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

//-- zung2r --------------------------------------------------------------------
void
LAPACK_DECL(zung2r)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zung2r");
    LAPACK_IMPL(zung2r)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zungbr");
    LAPACK_IMPL(zungbr)(VECT,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zunghr");
    LAPACK_IMPL(zunghr)(N,
                        ILO,
                        IHI,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

//-- zungl2 --------------------------------------------------------------------
void
LAPACK_DECL(zungl2)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zungl2");
    LAPACK_IMPL(zungl2)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zunglq");
    LAPACK_IMPL(zunglq)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zungql");
    LAPACK_IMPL(zungql)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zungqr");
    LAPACK_IMPL(zungqr)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

//-- zungr2 --------------------------------------------------------------------
void
LAPACK_DECL(zungr2)(const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zungr2");
    LAPACK_IMPL(zungr2)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zungrq");
    LAPACK_IMPL(zungrq)(M,
                        N,
                        K,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

//-- zungtr --------------------------------------------------------------------
void
LAPACK_DECL(zungtr)(const char               *UPLO,
                    const INTEGER            *N,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zungtr");
    LAPACK_IMPL(zungtr)(UPLO,
                        N,
                        A,
                        LDA,
                        TAU,
                        WORK,
                        LWORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zunm2l");
    LAPACK_IMPL(zunm2l)(SIDE,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zunm2r");
    LAPACK_IMPL(zunm2r)(SIDE,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zunmbr");
    LAPACK_IMPL(zunmbr)(VECT,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zunmhr");
    LAPACK_IMPL(zunmhr)(SIDE,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zunml2");
    LAPACK_IMPL(zunml2)(SIDE,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zunmlq");
    LAPACK_IMPL(zunmlq)(SIDE,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zunmql");
    LAPACK_IMPL(zunmql)(SIDE,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zunmqr");
    LAPACK_IMPL(zunmqr)(SIDE,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zunmr2");
    LAPACK_IMPL(zunmr2)(SIDE,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zunmr3");
    LAPACK_IMPL(zunmr3)(SIDE,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zunmrq");
    LAPACK_IMPL(zunmrq)(SIDE,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zunmrz");
    LAPACK_IMPL(zunmrz)(SIDE,
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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zunmtr");
    LAPACK_IMPL(zunmtr)(SIDE,
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

//-- zupgtr --------------------------------------------------------------------
void
LAPACK_DECL(zupgtr)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *AP,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *Q,
                    const INTEGER            *LDQ,
                    DOUBLE_COMPLEX           *WORK,
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zupgtr");
    LAPACK_IMPL(zupgtr)(UPLO,
                        N,
                        AP,
                        TAU,
                        Q,
                        LDQ,
                        WORK,
                        INFO);
}

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
                    INTEGER                  *INFO)
{
    DEBUG_LAPACK_STUB("zupmtr");
    LAPACK_IMPL(zupmtr)(SIDE,
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

//-- cgetrf --------------------------------------------------------------------
void
LAPACK_DECL(cgetrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    FLOAT_COMPLEX    *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("cgetrf");
    LAPACK_IMPL(cgetrf)(M,
                        N,
                        A,
                        LDA,
                        IPIV,
                        INFO);
}

//-- cgetrs --------------------------------------------------------------------
void
LAPACK_DECL(cgetrs)(const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    FLOAT_COMPLEX          *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    FLOAT_COMPLEX          *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("cgetrs");
    LAPACK_IMPL(cgetrs)(TRANS,
                        N,
                        NRHS,
                        A,
                        LDA,
                        IPIV,
                        B,
                        LDB,
                        INFO);
}

//-- clag2z --------------------------------------------------------------------
void
LAPACK_DECL(clag2z)(const INTEGER    *M,
                    const INTEGER    *N,
                    FLOAT_COMPLEX          *SA,
                    const INTEGER    *LDSA,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("clag2z");
    LAPACK_IMPL(clag2z)(M,
                        N,
                        SA,
                        LDSA,
                        A,
                        LDA,
                        INFO);
}

//-- cpotrf --------------------------------------------------------------------
void
LAPACK_DECL(cpotrf)(const char       *UPLO,
                    const INTEGER    *N,
                    FLOAT_COMPLEX          *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("cpotrf");
    LAPACK_IMPL(cpotrf)(UPLO,
                        N,
                        A,
                        LDA,
                        INFO);
}

//-- cpotrs --------------------------------------------------------------------
void
LAPACK_DECL(cpotrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    FLOAT_COMPLEX    *A,
                    const INTEGER    *LDA,
                    FLOAT_COMPLEX    *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("cpotrs");
    LAPACK_IMPL(cpotrs)(UPLO,
                        N,
                        NRHS,
                        A,
                        LDA,
                        B,
                        LDB,
                        INFO);
}

//-- ilaprec -------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilaprec)(const char   *PREC)
{
    DEBUG_LAPACK_STUB("ilaprec");
    return LAPACK_IMPL(ilaprec)(PREC);
}

//-- chla_transtype ------------------------------------------------------------
char
LAPACK_DECL(chla_transtype)(const INTEGER    *TRANS)
{
    DEBUG_LAPACK_STUB("chla_transtype");
    return LAPACK_IMPL(chla_transtype)(TRANS);
}

//-- claswp --------------------------------------------------------------------
void
LAPACK_DECL(claswp)(const INTEGER    *N,
                    FLOAT_COMPLEX    *A,
                    const INTEGER    *LDA,
                    const INTEGER    *K1,
                    const INTEGER    *K2,
                    const INTEGER    *IPIV,
                    const INTEGER    *INCX)
{
    DEBUG_LAPACK_STUB("claswp");
    LAPACK_IMPL(claswp)(N,
                        A,
                        LDA,
                        K1,
                        K2,
                        IPIV,
                        INCX);
}

//-- izmax1 --------------------------------------------------------------------
INTEGER
LAPACK_DECL(izmax1)(const INTEGER            *N,
                    const DOUBLE_COMPLEX     *CX,
                    const INTEGER            *INCX)
{
    DEBUG_LAPACK_STUB("izmax1");
    return LAPACK_IMPL(izmax1)(N,
                               CX,
                               INCX);
}

//-- ilazlc --------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilazlc)(const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA)
{
    DEBUG_LAPACK_STUB("ilazlc");
    return LAPACK_IMPL(ilazlc)(M,
                               N,
                               A,
                               LDA);
}

//-- ilazlr --------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilazlr)(const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA)
{
    DEBUG_LAPACK_STUB("ilazlr");
    return LAPACK_IMPL(ilazlr)(M,
                               N,
                               A,
                               LDA);
}

//-- cpotf2 --------------------------------------------------------------------
void
LAPACK_DECL(cpotf2)(const char       *UPLO,
                    const INTEGER    *N,
                    FLOAT_COMPLEX    *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    DEBUG_LAPACK_STUB("cpotf2");
    LAPACK_IMPL(cpotf2)(UPLO,
                        N,
                        A,
                        LDA,
                        INFO);
}

//-- clacgv --------------------------------------------------------------------
void
LAPACK_DECL(clacgv)(const INTEGER    *N,
                    FLOAT_COMPLEX    *X,
                    const INTEGER    *INCX)
{
    DEBUG_LAPACK_STUB("clacgv");
    LAPACK_IMPL(clacgv)(N,
                        X,
                        INCX);
}

