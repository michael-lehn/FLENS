#include <flens/lapack/interface/include/config.h>
#include <flens/flens.cxx>

namespace flens { namespace lapack {

extern "C" {

MachineParameter
getMachineParameter_(const char c)
{
    MachineParameter machineParameter = Eps;

    switch (c) {
        case 'E':
        case 'e': machineParameter = Eps;
                  break;
        case 'S':
        case 's': machineParameter = SafeMin;
                  break;
        case 'B':
        case 'b': machineParameter = Base;
                  break;
        case 'P':
        case 'p': machineParameter = Precision;
                  break;
        case 'N':
        case 'n': machineParameter = Mantissa;
                  break;
        case 'R':
        case 'r': machineParameter = Rounding;
                  break;
        case 'M':
        case 'm': machineParameter = UnderflowExp;
                  break;
        case 'U':
        case 'u': machineParameter = UnderflowThreshold;
                  break;
        case 'L':
        case 'l': machineParameter = OverflowExp;
                  break;
        case 'O':
        case 'o': machineParameter = OverflowThreshold;
                  break;
        default:  ASSERT(0);
    }
    return machineParameter;
}

//-- dlamch --------------------------------------------------------------------
// TODO: change this back to "LAPACK_DECL(dlamch)"
DOUBLE
LAPACK_DECL(__dlamch)(const char *c)
{
    return lamch<DOUBLE>(getMachineParameter_(*c));
}

} // extern "C"

} } // namespace lapack, flens
