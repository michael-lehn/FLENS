#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define LABAD         slabad_
#   define LABAD_REF     slabad
#   define LABAD_NAME    "SLABAD"
#elif DOUBLE
#   define LABAD         dlabad_
#   define LABAD_REF     dlabad
#   define LABAD_NAME    "DLABAD"
#elif COMPLEX_SINGLE
#   define LABAD         clabad_
#   define LABAD_REF     clabad
#   define LABAD_NAME    "CLABAD"
#elif COMPLEX_DOUBLE
#   define LABAD         zlabad_
#   define LABAD_REF     zlabad
#   define LABAD_NAME    "ZLABAD"
#endif

using namespace flens;

extern "C" {


void
LABAD_REF(FLOAT *SMALL, FLOAT *LARGE);


void
LABAD(FLOAT *SMALL, FLOAT *LARGE)
{
#   ifndef USE_LABAD_REF
    ASSERT(0);  // not yet implemented
#   else
    LABAD_REF(SMALL, LARGE);
#   endif
}

} // extern "C"

