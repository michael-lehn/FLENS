//
//  Regression test for lsame_ in the Fortran interface of FLENS-LAPACK
//  (flens/lapack/interface/src/lsame.cc, part of lapack_FLENS.a).
//
//  LSAME compares two characters regardless of case.  The FLENS version
//  compared case-sensitively.  As lapack_FLENS.a exports lsame_, Fortran
//  routines linked with it (e.g. the reference BLAS) rejected lower case
//  options: in the Netlib test suite (zgbak) ZGEMM('n', 'n', ...) reported
//  "parameter number 1 had an illegal value" and returned without computing
//  anything.
//

#include <flens/lapack/interface/src/lsame.cc>

#include <cxxstd/iostream.h>

using namespace std;

static char
upper(char c)
{
    return (c>='a' && c<='z') ? char(c-'a'+'A') : c;
}

int
main()
{
    int failures = 0;

    for (int a=32; a<127; ++a) {
        for (int b=32; b<127; ++b) {
            const char ca = char(a), cb = char(b);
            const bool expected = (upper(ca)==upper(cb));
            const bool result   = lsame_(&ca, &cb)!=0;
            if (result!=expected) {
                if (++failures<=10) {
                    cerr << "FAILED: lsame_('" << ca << "', '" << cb
                         << "') returned " << result << endl;
                }
            }
        }
    }

    if (failures) {
        cerr << failures << " check(s) FAILED" << endl;
        return 1;
    }
    cout << "lsame-interface: all checks passed" << endl;
    return 0;
}
