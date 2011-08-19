#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define LAPY2         slapy2_
#   define LAPY2_REF     slapy2
#   define LAPY2_NAME    "SLAPY2"
#elif DOUBLE
#   define LAPY2         dlapy2_
#   define LAPY2_REF     dlapy2
#   define LAPY2_NAME    "DLAPY2"
#elif COMPLEX_SINGLE
#   define LAPY2         clapy2_
#   define LAPY2_REF     clapy2
#   define LAPY2_NAME    "CLAPY2"
#elif COMPLEX_DOUBLE
#   define LAPY2         zlapy2_
#   define LAPY2_REF     zlapy2
#   define LAPY2_NAME    "ZLAPY2"
#endif

using namespace flens;

extern "C" {


FLOAT
LAPY2_REF(FLOAT *X, FLOAT*Y);


FLOAT
LAPY2(FLOAT *X, FLOAT*Y)
{
#   ifndef USE_LAPY2_REF

    double _x = *X, _y = *Y;

    //FLOAT flensResult = lapack::lapy2(*X, *Y);
    FLOAT flensResult = lapack::lapy2(_x, _y);
    FLOAT refResult = LAPY2_REF(X, Y);

    std::cerr << "flensResult = " << flensResult << std::endl;
    std::cerr << "refResult =   " << refResult << std::endl;

    FLOAT diff = refResult - flensResult;
    std::cerr << "diff =        " << diff << std::endl;

    ASSERT(flensResult==refResult);

    return flensResult;
#   else
    return LAPY2_REF(X, Y);
#   endif
}

} // extern "C"

