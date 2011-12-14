#include <flens/lapack/interface/include/config.h>

extern "C" {

//-- lsame ---------------------------------------------------------------------
LOGICAL
LAPACK_DECL(lsame)(const char       *CA,
                   const char       *CB)
{
    DEBUG_FLENS_LAPACK("lsame");
    return (*CA==*CB) ? LOGICAL(1) : LOGICAL(0);
}

} // extern "C"
