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
       SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
       SUBROUTINE ZLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LAEV2_TCC
#define FLENS_LAPACK_LA_LAEV2_TCC 1

#include <cxxstd/complex.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//
//  Real variant
//
template <typename T>
void
laev2_impl(const T &a,
           const T &b,
           const T &c,
           T       &rt1,
           T       &rt2,
           T       &cs1,
           T       &sn1)
{
    using cxxblas::pow;
    using std::abs;
    using std::sqrt;

    const T Zero(0), Half(0.5), One(1), Two(2);

//
//  Compute the eigenvalues
//
    T sm = a + c;
    T df = a - c;
    T adf = abs(df);
    T tb = b + b;
    T ab = abs(tb);

    T acmx, acmn;

    if (abs(a)>abs(c)) {
        acmx = a;
        acmn = c;
    } else {
        acmx = c;
        acmn = a;
    }

    T   rt;
    int sign1;

    if (adf>ab) {
        rt = adf*sqrt(One + pow(ab/adf, 2));
    } else if (adf<ab) {
        rt = ab*sqrt(One + pow(adf/ab, 2));
    } else {
//
//      Includes case AB=ADF=0
//
        rt = ab*sqrt(Two);
    }
    if (sm<Zero) {
        rt1   = Half*(sm-rt);
        sign1 = -1;
//
//      Order of execution important.
//      To get fully accurate smaller eigenvalue,
//      next line needs to be executed in higher precision.
//
        rt2 = (acmx/rt1)*acmn - (b/rt1)*b;
    } else if (sm>Zero) {
        rt1   = Half*(sm+rt);
        sign1 = 1;
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
        rt1   =  Half*rt;
        rt2   = -Half*rt;
        sign1 = 1;
    }
//
//  Compute the eigenvector
//
    T   cs;
    int sign2;

    if (df>=Zero) {
        cs = df + rt;
        sign2 = 1;
    } else {
        cs = df - rt;
        sign2 = -1;
    }

    T acs = abs(cs);

    if (acs>ab) {
        T ct  = -tb / cs;

        sn1 = One / sqrt(One + ct*ct);
        cs1 = ct*sn1;
    } else {
        if (ab==Zero) {
            cs1 = One;
            sn1 = Zero;
        } else {
            T tn = -cs / tb;

            cs1 = One / sqrt(One + tn*tn);
            sn1 = tn*cs1;
        }
    }
    if (sign1==sign2) {
        T tn = cs1;

        cs1 = -sn1;
        sn1 =  tn;
    }
}

//
//  Complex variant
//
template <typename T>
void
laev2_impl(const std::complex<T> &a,
           const std::complex<T> &b,
           const std::complex<T> &c,
           T                     &rt1,
           T                     &rt2,
           T                     &cs1,
           std::complex<T>       &sn1)
{
    using std::abs;
    using std::conj;
    using std::real;

    const T Zero(0), One(1);

    std::complex<T>  w;
    T                t;

    if (abs(b)==Zero) {
        w = One;
    } else {
        w = conj(b) / abs(b);
    }

    laev2(real(a), abs(b), real(c), rt1, rt2, cs1, t);
    sn1 = w*t;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//
//  Real variant
//
template <typename T>
void
laev2_impl(const T &a,
           const T &b,
           const T &c,
           T       &rt1,
           T       &rt2,
           T       &cs1,
           T       &sn1)
{
    cxxlapack::laev2(a, b, c, rt1, rt2, cs1, sn1);
}

//
//  Complex variant
//
template <typename T>
void
laev2_impl(const std::complex<T> &a,
           const std::complex<T> &b,
           const std::complex<T> &c,
           T                     &rt1,
           T                     &rt2,
           T                     &cs1,
           std::complex<T>       &sn1)
{
    cxxlapack::laev2(a, b, c, rt1, rt2, cs1, sn1);
}

} // namespace external

#endif

//== public interface ==========================================================

//
//  Real variant
//
template <typename T>
typename RestrictTo<IsReal<T>::value,
         void>::Type
laev2(const T &a,
      const T &b,
      const T &c,
      T       &rt1,
      T       &rt2,
      T       &cs1,
      T       &sn1)
{
    LAPACK_DEBUG_OUT("(d)laev2");
//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    T rt1_org  = rt1;
    T rt2_org  = rt2;
    T cs1_org  = cs1;
    T sn1_org  = sn1;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::laev2_impl(a, b, c, rt1, rt2, cs1, sn1);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    T rt1_generic  = rt1;
    T rt2_generic  = rt2;
    T cs1_generic  = cs1;
    T sn1_generic  = sn1;

    rt1 = rt1_org;
    rt2 = rt2_org;
    cs1 = cs1_org;
    sn1 = sn1_org;

    external::laev2_impl(a, b, c, rt1, rt2, cs1, sn1);

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
    if (! isIdentical(cs1_generic, cs1, "cs1_generic", "cs1")) {
        std::cerr << "CXXLAPACK: cs1_generic = " << cs1_generic << std::endl;
        std::cerr << "F77LAPACK: cs1 = " << cs1 << std::endl;
        failed = true;
    }
    if (! isIdentical(sn1_generic, sn1, "sn1_generic", "sn1")) {
        std::cerr << "CXXLAPACK: sn1_generic = " << sn1_generic << std::endl;
        std::cerr << "F77LAPACK: sn1 = " << sn1 << std::endl;
        failed = true;
    }
    if (failed) {
        ASSERT(0);
    }
#   endif
}

//
//  Complex variant
//
template <typename T>
typename RestrictTo<IsReal<T>::value,
         void>::Type
laev2(const std::complex<T> &a,
      const std::complex<T> &b,
      const std::complex<T> &c,
      T                     &rt1,
      T                     &rt2,
      T                     &cs1,
      std::complex<T>       &sn1)
{
    LAPACK_DEBUG_OUT("(z)laev2");
//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    T               rt1_org  = rt1;
    T               rt2_org  = rt2;
    T               cs1_org  = cs1;
    std::complex<T> sn1_org  = sn1;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::laev2_impl(a, b, c, rt1, rt2, cs1, sn1);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    T               rt1_generic  = rt1;
    T               rt2_generic  = rt2;
    T               cs1_generic  = cs1;
    std::complex<T> sn1_generic  = sn1;

    rt1 = rt1_org;
    rt2 = rt2_org;
    cs1 = cs1_org;
    sn1 = sn1_org;

    external::laev2_impl(a, b, c, rt1, rt2, cs1, sn1);

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
    if (! isIdentical(cs1_generic, cs1, "cs1_generic", "cs1")) {
        std::cerr << "CXXLAPACK: cs1_generic = " << cs1_generic << std::endl;
        std::cerr << "F77LAPACK: cs1 = " << cs1 << std::endl;
        failed = true;
    }
    if (! isIdentical(sn1_generic, sn1, "sn1_generic", "sn1")) {
        std::cerr << "CXXLAPACK: sn1_generic = " << sn1_generic << std::endl;
        std::cerr << "F77LAPACK: sn1 = " << sn1 << std::endl;
        failed = true;
    }
    if (failed) {
        ASSERT(0);
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAEV2_TCC
