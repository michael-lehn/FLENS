#include <cxxstd/cassert.h>
#include <cxxstd/cctype.h>

#include <flens/lapack/interface/include/config.h>

#ifdef LAPACK_IMPL
#   undef   LAPACK_IMPL
#endif
#define  LAPACK_IMPL(x)     x

extern "C" {

INTEGER
LAPACK_DECL(ilaenv)(const INTEGER *SPEC,
                    const char *NAME,
                    const char *OPTS,
                    const INTEGER *N1,
                    const INTEGER *N2,
                    const INTEGER *N3,
                    const INTEGER *N4,
                    int NAME_LEN,
                    int OPTS_LEN);

//
// for "ilaenv" and "xerbla" LAPACK_IMPL(..) calls the LAPACK_DECL(..)
// such that the routines from the LAPACK test suite get called
//
INTEGER
LAPACK_IMPL(ilaenv)(const INTEGER *SPEC,
                    const char *NAME,
                    const char *OPTS,
                    const INTEGER *N1,
                    const INTEGER *N2,
                    const INTEGER *N3,
                    const INTEGER *N4,
                    int NAME_LEN,
                    int OPTS_LEN)
{
    return LAPACK_DECL(ilaenv)(SPEC,
                               NAME,
                               OPTS,
                               N1,
                               N2,
                               N3,
                               N4,
                               NAME_LEN,
                               OPTS_LEN);
}

void
LAPACK_IMPL(xerbla)(const char      *SRNAME,
                    const INTEGER   *INFO,
                    int             SRNAME_LEN)
{
    return LAPACK_DECL(xerbla)(SRNAME,
                               INFO,
                               SRNAME_LEN);
}

} // extern "C"
