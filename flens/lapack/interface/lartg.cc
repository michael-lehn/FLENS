#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define LARTG         slartg_
#   define LARTG_REF     slartg
#   define LARTG_NAME    "SLARTG"
#elif DOUBLE
#   define LARTG         dlartg_
#   define LARTG_REF     dlartg
#   define LARTG_NAME    "DLARTG"
#elif COMPLEX_SINGLE
#   define LARTG         clartg_
#   define LARTG_REF     clartg
#   define LARTG_NAME    "CLARTG"
#elif COMPLEX_DOUBLE
#   define LARTG         zlartg_
#   define LARTG_REF     zlartg
#   define LARTG_NAME    "ZLARTG"
#endif

using namespace flens;

extern "C" {


void
LARTG_REF(FLOAT *F, FLOAT *G, FLOAT *CS, FLOAT *SN, FLOAT *R);


void
LARTG(FLOAT *F, FLOAT *G, FLOAT *CS, FLOAT *SN, FLOAT *R)
{
#   ifndef USE_LARTG_REF
    ASSERT(0);  // not yet implemented
#   else
    LARTG_REF(F, G, CS, SN, R);
#   endif
}

} // extern "C"

