#include <flens/lapack/interface/interface.h>


#ifdef SINGLE
#   define LARNV         slarnv_
#   define LARNV_REF     slarnv
#   define LARNV_NAME    "SLARNV"
#elif DOUBLE
#   define LARNV         dlarnv_
#   define LARNV_REF     dlarnv
#   define LARNV_NAME    "DLARNV"
#elif COMPLEX_SINGLE
#   define LARNV         clarnv_
#   define LARNV_REF     clarnv
#   define LARNV_NAME    "CLARNV"
#elif COMPLEX_DOUBLE
#   define LARNV         zlarnv_
#   define LARNV_REF     zlarnv
#   define LARNV_NAME    "ZLARNV"
#endif


extern "C" {


void
LARNV_REF(const INT *IDIST, INT *ISEED, const INT *N, FLOAT *X);


void
LARNV(const INT *IDIST, INT *ISEED, const INT *N, FLOAT *X)
{
#   ifndef USE_LARNV_REF
    ASSERT(0);  // not yet implemented
#   else
    LARNV_REF(IDIST, ISEED, N, X);
#   endif
}

} // extern "C"

