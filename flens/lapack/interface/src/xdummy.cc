#include <cassert>
#include <cctype>

#include <flens/lapack/interface/include/config.h>

#define ASSERT(x)   assert(x)

extern "C" {

#include <flens/lapack/interface/include/dummy.in.cc>


//
// for "ilaenv" we have to reverse the pass-through such that the
// ilaenv version from the testsuit gets called
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
