/*
 *   Copyright (c) 2011, Michael Lehn
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
       SUBROUTINE DLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
 *
 *  -- LAPACK auxiliary routine (version 3.2.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     June 2010
 */

#ifndef FLENS_LAPACK_EIG_LANV2_TCC
#define FLENS_LAPACK_EIG_LANV2_TCC 1

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename T>
void
lanv2_generic(T &a, T &b, T &c, T &d,
              T &rt1r, T &rt1i,
              T &rt2r, T &rt2i,
              T &cs, T &sn)
{
    using std::abs;
    using std::max;
    using std::min;
    using std::swap;

    const T Zero(0), Half(0.5), One(1), Multpl(4);

    const T eps = lamch<T>(Precision);
    if (c==Zero) {
        cs = One;
        sn = Zero;
    } else if (b==Zero) {
//
//      Swap rows and columns
//
        cs = Zero;
        sn = One;
        swap(a, d);
        b = -c;
        c = Zero;
    } else if ((a-d==Zero) && (sign(One,b)!=sign(One,c))) {
        cs = One;
        sn = Zero;
    } else {
        T temp = a-d;
        const T p = Half*temp;
        const T bcMax = max(abs(b), abs(c));
        const T bcMis = min(abs(b), abs(c)) * sign(One,b) * sign (One, c);
        const T scale = max(abs(p), bcMax);
        T z = (p/scale)*p + (bcMax/scale)*bcMis;
//
//      If Z is of the order of the machine accuracy, postpone the
//      decision on the nature of eigenvalues
//
        if (z>=Multpl*eps) {
//
//          Real eigenvalues. Compute A and D.
//
            z = p + sign(sqrt(scale)*sqrt(z), p);
            a = d + z;
            d = d - (bcMax/z)*bcMis;
//
//          Compute B and the rotation matrix
//
            const T tau = lapy2(c, z);
            cs = z/tau;
            sn = c/tau;
            b = b -c;
            c = Zero;
        } else {
//
//          Complex eigenvalues, or real (almost) equal eigenvalues.
//          Make diagonal elements equal.
//
            const T sigma = b + c;
            const T tau = lapy2(sigma, temp);
            cs = sqrt(Half*(One+abs(sigma)/tau));
            sn = -(p/(tau*cs)) * sign(One, sigma);
//
//          Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
//                  [ CC  DD ]   [ C  D ] [ SN  CS ]
//
            const T aa =  a*cs + b*sn;
            const T bb = -a*sn + b*cs;
            const T cc =  c*cs + d*sn;
            const T dd = -c*sn + d*cs;
//
//          Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
//                  [ C  D ]   [-SN  CS ] [ CC  DD ]
//
            a =  aa*cs + cc*sn;
            b =  bb*cs + dd*sn;
            c = -aa*sn + cc*cs;
            d = -bb*sn + dd*cs;

            temp = Half * (a+d);
            a = temp;
            d = temp;

            if (c!=Zero) {
                if (b!=Zero) {
                    if (sign(One, b)==sign(One,c)) {
//
//                      Real eigenvalues: reduce to upper triangular form
//
                        const T sab = sqrt(abs(b));
                        const T sac = sqrt(abs(c));
                        const T p = sign(sab*sac, c);
                        const T tau = One / sqrt(abs(b+c));
                        a = temp + p;
                        d = temp - p;
                        b = b -c;
                        c = Zero;
                        const T cs1 = sab*tau;
                        const T sn1 = sac*tau;
                        temp = cs*cs1 - sn*sn1;
                        sn = cs*sn1 + sn*cs1;
                        cs = temp;
                    }
               } else {
                    b = -c;
                    c = Zero;
                    temp = cs;
                    cs = -sn;
                    sn = temp;
               }
            }
        }
    }
//
//  Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
//
    rt1r = a;
    rt2r = d;
    if (c==Zero) {
        rt1i = Zero;
        rt2i = Zero;
    } else {
        rt1i = sqrt(abs(b))*sqrt(abs(c));
        rt2i = -rt1i;
    }
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename T>
void
lanv2_native(T &a, T &b, T &c, T &d,
             T &rt1r, T &rt1i,
             T &rt2r, T &rt2i,
             T &cs, T &sn)
{
    if (IsSame<T,DOUBLE>::value) {
        LAPACK_IMPL(dlanv2)(&a, &b, &c, &d,
                            &rt1r, &rt1i,
                            &rt2r, &rt2i,
                            &cs, &sn);
    } else {
        ASSERT(0);
    }
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename T>
void
lanv2(T &a, T &b, T &c, T &d, T &rt1r, T &rt1i, T &rt2r, T &rt2i, T &cs, T &sn)
{
    LAPACK_DEBUG_OUT("lanv2");

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    T _a    = a;
    T _b    = b;
    T _c    = c;
    T _d    = d;
    T _rt1r = rt1r;
    T _rt1i = rt1i;
    T _rt2r = rt2r;
    T _rt2i = rt2i;
    T _cs   = cs;
    T _sn   = sn;
#   endif

//
//  Call implementation
//
    lanv2_generic(a, b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    lanv2_native(_a, _b, _c, _d, _rt1r, _rt1i, _rt2r, _rt2i, _cs, _sn);

    bool failed = false;
    if (! isIdentical(a, _a, " a", "_a")) {
        std::cerr << "CXXLAPACK:  a = " << a << std::endl;
        std::cerr << "F77LAPACK: _a = " << _a << std::endl;
        failed = true;
    }

    if (! isIdentical(b, _b, " b", "_b")) {
        std::cerr << "CXXLAPACK:  b = " << b << std::endl;
        std::cerr << "F77LAPACK: _b = " << _b << std::endl;
        failed = true;
    }

    if (! isIdentical(c, _c, " c", "_c")) {
        std::cerr << "CXXLAPACK:  c = " << c << std::endl;
        std::cerr << "F77LAPACK: _c = " << _c << std::endl;
        failed = true;
    }

    if (! isIdentical(d, _d, " d", "_d")) {
        std::cerr << "CXXLAPACK:  d = " << d << std::endl;
        std::cerr << "F77LAPACK: _d = " << _d << std::endl;
        failed = true;
    }

    if (! isIdentical(rt1r, _rt1r, " rt1r", "_rt1r")) {
        std::cerr << "CXXLAPACK:  rt1r = " << rt1r << std::endl;
        std::cerr << "F77LAPACK: _rt1r = " << _rt1r << std::endl;
        failed = true;
    }

    if (! isIdentical(rt1i, _rt1i, " rt1i", "_rt1i")) {
        std::cerr << "CXXLAPACK:  rt1i = " << rt1i << std::endl;
        std::cerr << "F77LAPACK: _rt1i = " << _rt1i << std::endl;
        failed = true;
    }

    if (! isIdentical(rt2r, _rt2r, " rt2r", "_rt2r")) {
        std::cerr << "CXXLAPACK:  rt2r = " << rt2r << std::endl;
        std::cerr << "F77LAPACK: _rt2r = " << _rt2r << std::endl;
        failed = true;
    }

    if (! isIdentical(rt2i, _rt2i, " rt2i", "_rt2i")) {
        std::cerr << "CXXLAPACK:  rt2i = " << rt2i << std::endl;
        std::cerr << "F77LAPACK: _rt2i = " << _rt2i << std::endl;
        failed = true;
    }

    if (! isIdentical(cs, _cs, " cs", "_cs")) {
        std::cerr << "CXXLAPACK:  cs = " << cs << std::endl;
        std::cerr << "F77LAPACK: _cs = " << _cs << std::endl;
        failed = true;
    }

    if (! isIdentical(sn, _sn, " sn", "_sn")) {
        std::cerr << "CXXLAPACK:  sn = " << sn << std::endl;
        std::cerr << "F77LAPACK: _sn = " << _sn << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_LANV2_TCC
