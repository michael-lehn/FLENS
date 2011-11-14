#include <flens/lapack/interface/include/config.h>

extern "C" {

//-- lsame ---------------------------------------------------------------------
LOGICAL
LAPACK_DECL(lsame)(const char       *CA,
                   const char       *CB)
{
    return (*CA==*CB) ? LOGICAL(1) : LOGICAL(0);
}

} // extern "C"