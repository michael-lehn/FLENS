/*
 *   Copyright (c) 2014, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* Based on
 *
       SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LAE2_TCC
#define FLENS_LAPACK_LA_LAE2_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename T>
void
lae2_impl(const T &a,
          const T &b,
          const T &c,
          T       &rt1,
          T       &rt2)
{
    using std::abs;
    using std::sqrt;

    const T Zero(0), Half(0.5), One(1), Two(2);

    T ab, acmn, acmx, adf, df, rt, sm, tb;

//
//  Compute the eigenvalues
//
    sm = a + c;
    df = a - c;
    adf = abs(df);
    tb = b + b;
    ab = abs(tb);

    if (abs(a)>abs(c)) {
        acmx = a;
        acmn = c;
    } else {
        acmx = c;
        acmn = a;
    }
    if (adf>ab) {
        rt = adf*sqrt(One+pow(ab/adf, 2));
    } else if (adf<ab) {
        rt = ab*sqrt(One+pow(adf/ab, 2));
    } else {
//
//      Includes case AB=ADF=0
//
        rt = ab*sqrt(Two);
    }
    if (sm<Zero) {
        rt1 = Half*(sm-rt);
//
//      Order of execution important.
//      To get fully accurate smaller eigenvalue,
//      next line needs to be executed in higher precision.
//
        rt2 = (acmx/rt1)*acmn - (b/rt1)*b;
    } else if (sm>Zero) {
        rt1 = Half*(sm+rt);
//
//      Order of execution important.
//      To get fully accurate smaller eigenvalue,
//      next line needs to be executed in higher precision.
//
        rt2 = (acmx/rt1)*acmn - (b/rt1)*b;
    } else {
//
//      Includes case RT1 = RT2 = 0
//
        rt1 =  Half*rt;
        rt2 = -Half*rt;
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename T>
void
lae2_impl(const T &a,
          const T &b,
          const T &c,
          T       &rt1,
          T       &rt2)
{
    cxxlapack::lae2(a, b, c, rt1, rt2);
}

} // namespace external

#endif

//== public interface ==========================================================

template <typename T>
void
lae2(const T &a,
     const T &b,
     const T &c,
     T       &rt1,
     T       &rt2)
{
//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    T rt1_org  = rt1;
    T rt2_org  = rt2;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::lae2_impl(a, b, c, rt1, rt2);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    T rt1_generic  = rt1;
    T rt2_generic  = rt2;

    rt1 = rt1_org;
    rt2 = rt2_org;

    external::lae2_impl(a, b, c, rt1, rt2);

    bool failed = false;
    if (! isIdentical(rt1_generic, rt1, "rt1_generic", "rt1")) {
        std::cerr << "CXXLAPACK: rt1_generic = " << rt1_generic << std::endl;
        std::cerr << "F77LAPACK: rt1 = " << rt1 << std::endl;
        failed = true;
    }
    if (! isIdentical(rt2_generic, rt2, "rt2_generic", "rt2")) {
        std::cerr << "CXXLAPACK: rt2_generic = " << rt2_generic << std::endl;
        std::cerr << "F77LAPACK: rt2 = " << rt2 << std::endl;
        failed = true;
    }
    if (failed) {
        ASSERT(0);
    }
#   endif

}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAE2_TCC
