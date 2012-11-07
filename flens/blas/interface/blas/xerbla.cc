#include <flens/blas/interface/blas/config.h>

extern "C" {

void
BLAS(xerbla)(const char *SRNAME, const INTEGER *INFO)
{
    std::cout << " ** On entry to " << SRNAME << " parameter number ";
    std::cout << *INFO << " had an illegal value" << std::endl;
}

} // extern "C
