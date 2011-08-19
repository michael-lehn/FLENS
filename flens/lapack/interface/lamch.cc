#include <flens/lapack/interface/interface.h>
#include <cassert>
#include <cctype>
#include <cstdio>

#ifdef SINGLE
#   define LAMCH         slamch_
#   define LAMCH_REF     slamch
#   define LAMCH_NAME    "SLAMCH"
#elif DOUBLE
#   define LAMCH         dlamch_
#   define LAMCH_REF     dlamch
#   define LAMCH_NAME    "DLAMCH"
#elif COMPLEX_SINGLE
#   define LAMCH         clamch_
#   define LAMCH_REF     clamch
#   define LAMCH_NAME    "CLAMCH"
#elif COMPLEX_DOUBLE
#   define LAMCH         zlamch_
#   define LAMCH_REF     zlamch
#   define LAMCH_NAME    "ZLAMCH"
#endif

using namespace flens;

extern "C" {

FLOAT
LAMCH_REF(const char *c);


FLOAT
LAMCH(const char *c)
{
    /*
        Relative machine underflow is taken to be    0.117549E-37
        Relative machine overflow  is taken to be    0.340282E+39
        Relative machine precision is taken to be    0.596046E-07
    */

#   ifndef USE_LAMCH_REF
    lapack::MachineParameter machineParameter;
    switch (*c) {
        case 'E':
        case 'e': machineParameter = lapack::Eps;
                  break;
        case 'S':
        case 's': machineParameter = lapack::SafeMin;
                  break;
        case 'B':
        case 'b': machineParameter = lapack::Base;
                  break;
        case 'P':
        case 'p': machineParameter = lapack::Precision;
                  break;
        case 'N':
        case 'n': machineParameter = lapack::Mantissa;
                  break;
        case 'R':
        case 'r': machineParameter = lapack::Rounding;
                  break;
        case 'M':
        case 'm': machineParameter = lapack::UnderflowExp;
                  break;
        case 'U':
        case 'u': machineParameter = lapack::UnderflowThreshold;
                  break;
        case 'L':
        case 'l': machineParameter = lapack::OverflowExp;
                  break;
        case 'O':
        case 'o': machineParameter = lapack::OverflowThreshold;
                  break;
        default:
            assert(0);
    }
    return lapack::lamch<FLOAT>(machineParameter);
#   else
    return LAMCH_REF(c);
#endif
}

} // extern "C"