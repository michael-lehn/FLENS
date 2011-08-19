#include <iostream>

#ifdef SINGLE
#   define FLOAT         float
#   define LAPY2         slapy2_
#   define LAPY2_REF     slapy2
#   define LAPY2_NAME    "SLAPY2"
#elif DOUBLE
#   define FLOAT         double
#   define LAPY2         dlapy2_
#   define LAPY2_REF     dlapy2
#   define LAPY2_NAME    "DLAPY2"
#elif COMPLEX_SINGLE
#   define FLOAT         float
#   define LAPY2         clapy2_
#   define LAPY2_REF     clapy2
#   define LAPY2_NAME    "CLAPY2"
#elif COMPLEX_DOUBLE
#   define FLOAT         double
#   define LAPY2         zlapy2_
#   define LAPY2_REF     zlapy2
#   define LAPY2_NAME    "ZLAPY2"
#endif

extern "C" {

FLOAT
LAPY2(FLOAT *X, FLOAT*Y);

}

int
main()
{
        
    FLOAT x = 2.26040777926887787030452159342530649155378341674804687500e-01;
    FLOAT y = 9.25487851391693094882384684751741588115692138671875000e-02;

    std::cerr << std::scientific;

    std::cerr.precision(80);
    std::cerr << "x = " << x << std::endl;
    FLOAT res = LAPY2(&x, &y);
    std::cerr << "LAPY2(&x, &y) = " << res << std::endl;
}