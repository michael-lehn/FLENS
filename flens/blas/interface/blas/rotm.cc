#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(srotm)(const INTEGER   *N,
            float           *X,
            const INTEGER   *INCX,
            float           *Y,
            const INTEGER   *INCY,
            const float     *P)
{
    ASSERT(0);
}

void
BLAS(drotm)(const INTEGER   *N,
            double          *X,
            const INTEGER   *INCX,
            double          *Y,
            const INTEGER   *INCY,
            const double    *P)
{
    ASSERT(0);
}

} // extern "C"
