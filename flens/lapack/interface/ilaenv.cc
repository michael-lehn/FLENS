#include <flens/lapack/interface/interface.h>

#define ILAENV_REF      ilaenv


extern "C" {

int
ILAENV_REF(int *SPEC, const char *NAME, const char *OPTS,
           int *N1, int *N2, int *N3, int *N4)
{
    // call the ilaenv function from the lapack test suite
    return ILAENV(SPEC, NAME, OPTS, N1, N2, N3, N4);
}

} // extern "C"

