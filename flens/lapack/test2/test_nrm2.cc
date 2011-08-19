#include <iostream>
#include <cmath>

#define INT     int

#ifdef SINGLE
#   define FLOAT         float
#   define NRM2          snrm2_
#   define NRM2_REF      snrm2
#   define NRM2_NAME     "SNRM2"
#elif DOUBLE
#   define FLOAT         double
#   define NRM2          dnrm2_
#   define NRM2_REF      dnrm2
#   define NRM2_NAME     "DNRM2"
#elif COMPLEX_SINGLE
#   define FLOAT         float
#   define NRM2          cnrm2_
#   define NRM2_REF      cnrm2
#   define NRM2_NAME     "CNRM2"
#elif COMPLEX_DOUBLE
#   define FLOAT         double
#   define NRM2          znrm2_
#   define NRM2_REF      znrm2
#   define NRM2_NAME     "ZNRM2"
#endif

extern "C" {

FLOAT
NRM2_REF(INT *N, FLOAT *X, INT *INCX);

}

template <typename IndexType, typename X, typename T>
void
nrm2_generic(IndexType n, const X *x, IndexType incX, T &norm)
{
    if (n<1) {
        norm = 0;
    } else if (n==1) {
        norm = abs(*x);
    } else {
        T scale = 0;
        T ssq = 1;
//      The following loop is equivalent to this call to the LAPACK
//      auxiliary routine:
//      CALL DLASSQ( N, X, INCX, SCALE, SSQ )
//
        for (IndexType i=0, iX=0; i<n; ++i, iX+=incX) {
            if (x[iX]!=T(0)) {
                T absXi = std::abs(x[iX]);
                if (scale<absXi) {
                    // TODO: use a sqr function
                    //ssq = T(1) + ssq * (scale/absXi) * (scale/absXi);
                    ssq = T(1) + ssq * pow(scale/absXi, 2);
                    scale = absXi;
                } else {
                    // TODO: use a sqr function
                    //ssq += (absXi/scale) * (absXi/scale);
                    ssq += pow(absXi/scale, 2);
                }
            }
        }
        norm = scale*sqrt(ssq);
    }
}



int
main()
{
    INT n = 9;
    INT incX = 1;
    FLOAT x[9];

    x[0] = -9.027226160473944549522684610565192997455596923828125000000000e-02;
    x[1] = -1.546738812154063724690900016867090016603469848632812500000000e-02;
    x[2] = -1.196703376118133932415332765231141820549964904785156250000000e-01;
    x[3] = -9.844044886088371137944363908900413662195205688476562500000000e-02;
    x[4] = -1.814500417583024549639247879895265214145183563232421875000000e-02;
    x[5] =  6.937052244473329576379683203413151204586029052734375000000000e-02;
    x[6] =  6.946503314396984540834978361090179532766342163085937500000000e-02;
    x[7] = -8.916166830333621673432276111270766705274581909179687500000000e-02;
    x[8] = -2.789078367624207854813533913329592905938625335693359375000000e-02;


    std::cerr << std::scientific;

    std::cerr.precision(80);
    FLOAT res = NRM2_REF(&n, x, &incX);
    std::cerr << "NRM2(&x, &y) = " << res << std::endl;
    
    double flensRes;
    nrm2_generic(n, x, incX, flensRes);

    std::cerr << "CXXNRM2 =      " << flensRes << std::endl;

    long double flensRes2;
    nrm2_generic(n, x, incX, flensRes2);

    std::cerr << "CXXNRM2 =      " << flensRes2 << std::endl;
    
}