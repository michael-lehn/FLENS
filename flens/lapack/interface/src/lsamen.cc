#include <flens/lapack/interface/include/config.h>

extern "C" {

//-- lsamen --------------------------------------------------------------------
LOGICAL
LAPACK_DECL(lsamen)(const INTEGER    *N,
                    const char       *CA,
                    const char       *CB)
{
    DEBUG_FLENS_LAPACK("lsamen");
    for (int i=0; i<*N; ++i) {
        if (toupper(CA[i])!=toupper(CB[i])) {
            return LOGICAL(0);
        }
    }
    return LOGICAL(1);
}

} // extern "C"
