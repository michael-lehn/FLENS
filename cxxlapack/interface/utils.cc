#include <cassert>
#include <cctype>
#include <cstdio>
#include <cxxlapack/cxxlapack.cxx>
#include <cxxblas/cxxblas.cxx>


#define LSAMEN  lsamen_
#define LSAME   lsame_
#define ILAVER  ilaver_
#define SECOND  second_

#ifdef SINGLE
#   define FLOAT   float
#   define LAMCH   slamch_
#endif


#ifdef SINGLE
#   define LAMCH_REF   slamch
#endif

extern "C" {

using cxxlapack::lamch;

void
ILAVER(int *major, int *minor, int *patch)
{
    *major = 3;
    *minor = 0;
    *patch = 0;
}

int
LSAMEN(const int *n, const char *str1, const char *str2)
{
    for (int i=0; i<*n; ++i) {
        if (toupper(str1[i])!=toupper(str2[i])) {
            return 0;
        }
    }
    return 1;
}

int
LSAME(const char *c1, const char *c2)
{
    return (toupper(*c1)==toupper(*c2)) ? 1 : 0;
}

float
SECOND()
{
    return 0;
}


#ifdef USE_LAMCH_REF
FLOAT
LAMCH_REF(const char *c);
#endif

FLOAT
LAMCH(const char *c)
{
    /*
        Relative machine underflow is taken to be    0.117549E-37
        Relative machine overflow  is taken to be    0.340282E+39
        Relative machine precision is taken to be    0.596046E-07
    */

#   ifdef USE_LAMCH_REF
    FLOAT value_ = LAMCH_REF(c);
#   endif

    FLOAT value = 0;
    switch (*c) {
        case 'E':
        case 'e': value = lamch<FLOAT>(cxxlapack::Eps);
                  break;
        case 'S':
        case 's': value = lamch<FLOAT>(cxxlapack::SafeMin);
                  break;
        case 'B':
        case 'b': value = lamch<FLOAT>(cxxlapack::Base);
                  break;
        case 'P':
        case 'p': value = lamch<FLOAT>(cxxlapack::Precision);
                  break;
        case 'N':
        case 'n': value = lamch<FLOAT>(cxxlapack::Mantissa);
                  break;
        case 'R':
        case 'r': value = lamch<FLOAT>(cxxlapack::Rounding);
                  break;
        case 'M':
        case 'm': value = lamch<FLOAT>(cxxlapack::UnderflowExp);
                  break;
        case 'U':
        case 'u': value = lamch<FLOAT>(cxxlapack::UnderflowThreshold);
                  break;
        case 'L':
        case 'l': value = lamch<FLOAT>(cxxlapack::OverflowExp);
                  break;
        case 'O':
        case 'o': value = lamch<FLOAT>(cxxlapack::OverflowThreshold);
                  break;
        default:
            assert(0);
    }

#   ifdef USE_LAMCH_REF
    if (value!=value_) {
        fprintf(stderr, "LAMCH:     %c -> %e\n", *c, value);
        fprintf(stderr, "LAMCH_REF: %c -> %e\n", *c, value_);
    }
#   endif

    return value;
}

} // extern "C"