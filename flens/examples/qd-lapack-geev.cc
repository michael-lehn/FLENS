#include <iostream>

///
///  Include the qd-headers
///
#include <qd/qd_real.h>
#include <qd/fpu.h>

#include <flens/flens.cxx>

using namespace std;
using namespace flens;

///
///  For using quad-double precision in this example change the typedef to
///
typedef qd_real   T;

///
///  For understanding the next code lines we first take a look at the
///  __QD Library__ documentation:
///
///  *--[BOX]------------------------------------------------------------*
///  |                                                                   |
///  | The algorithms in the QD library assume IEEE double precision     |
///  | floating point arithmetic. Since Intel x86 processors have        |
///  | extended (80-bit) floating point registers, the round-to-double   |
///  | flag must be enabled in the control word of the FPU for this      |
///  | library to function properly under x86 processors. The function   |
///  | `fpu_fix_start` turns on the round-to-double bit in the FPU       |
///  | control word, while `fpu_fix_end` will restore the original       |
///  | state.                                                            |
///  |                                                                   |
///  *-------------------------------------------------------------------*
///
///  So the first thing we do in main is turning on the correct rounding ...
///
int
main()
{
    unsigned int old_cw;
    fpu_fix_start(&old_cw);

    typedef GeMatrix<FullStorage<T, ColMajor> >   Matrix;
    typedef DenseVector<Array<T> >                Vector;

    const int n = 5;

    Matrix   A(n, n), VL(n, n), VR(n, n);
    Vector   wr(n), wi(n);


    A =  2,   3,  -1,   0,  2,
        -6,  -5,   0,   2, -6,
         2,  -5,   6,  -6,  2,
         2,   3,  -1,   0,  8,
        -6,  -5,  10,   2, -6;

    cerr << "A = " << A << endl;

    Vector   work;
    lapack::ev(true, true, A, wr, wi, VL, VR, work);

    cerr << "wr = " << wr << endl;
    cerr << "wi = " << wi << endl;
    cerr << "VL = " << VL << endl;
    cerr << "VR = " << VR << endl;

    ///
    ///  ... and at the end restore FPU rounding behavior as mentioned above.
    ///
    fpu_fix_end(&old_cw);
}

///
///  :links: __QD Library__ -> http://crd-legacy.lbl.gov/~dhbailey/mpdist/
///
