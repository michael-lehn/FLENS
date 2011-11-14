//-- xerbla --------------------------------------------------------------------
void
LAPACK_DECL(xerbla)(const char      *SRNAME,
                    const INTEGER   *INFO,
                    int             SRNAME_LEN);

//-- ilaenv --------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilaenv)(const INTEGER   *SPEC,
                    const char      *NAME,
                    const char      *OPTS,
                    const INTEGER   *N1,
                    const INTEGER   *N2,
                    const INTEGER   *N3,
                    const INTEGER   *N4,
                    int             NAME_LEN,
                    int             OPTS_LEN);
