#include <cxxstd/locale.h>
#include <flens/lapack/interface/include/config.h>

extern "C" {

//-- lsamen --------------------------------------------------------------------
LOGICAL
LAPACK_DECL(lsamen)(const INTEGER    *N,
                    const char       *CA,
                    const char       *CB)
{
    std::locale loc;
    for (int i=0; i<*N; ++i) {
        if (std::toupper(CA[i], loc)!=std::toupper(CB[i],loc)) {
            return LOGICAL(0);
        }
    }
    return LOGICAL(1);
}

} // extern "C"
