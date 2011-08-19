#include <flens/lapack/interface/interface.h>

#define DSECND          dsecnd_
#define DSECND_REF      dsecnd
#define DSECND_NAME     "DSECND"

using namespace flens;

extern "C" {


double
DSECND_REF();


double
DSECND()
{
#   ifndef USE_LABAD_REF
    ASSERT(0);  // not yet implemented
#   else
    return DSECND_REF();
#   endif
}

} // extern "C"

