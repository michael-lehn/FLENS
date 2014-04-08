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
       SUBROUTINE DLARTG( F, G, CS, SN, R )
       SUBROUTINE ZLARTG( F, G, CS, SN, R )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LARTG_TCC
#define FLENS_LAPACK_LA_LARTG_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//
//  Real variant
//
template <typename T>
typename RestrictTo<IsReal<T>::value,
         void>::Type
lartg_impl(const T &f, const T &g, T &cs, T &sn, T &r)
{
    using std::abs;
    using std::max;
    using cxxblas::pow;

    const T     Zero(0), One(1), Two(2);
    static bool first = true;

    static T   eps, safeMin, safeMin2, safeMax2;

    if (first) {
        safeMin     = lamch<T>(SafeMin);
        eps         = lamch<T>(Eps);

        const int exp = explicit_cast<T,int>(log(safeMin/eps)
                                             / log(lamch<T>(Base))
                                             / Two);
        safeMin2    = pow(lamch<T>(Base), exp);
        safeMax2    = One / safeMin2;
        first = false;
    }

    if (g==Zero) {
        cs = One;
        sn = Zero;
        r = f;
    } else if (f==Zero) {
        cs = Zero;
        sn = One;
        r = g;
    } else {
        T f1 = f;
        T g1 = g;
        T scale = max(abs(f1), abs(g1));
        if (scale>=safeMax2) {
            int count = 0;
            do {
                ++count;
                f1 *= safeMin2;
                g1 *= safeMin2;
                scale = max(abs(f1), abs(g1));
            } while (scale>=safeMax2);
            r = sqrt(pow(f1,2) + pow(g1,2));
            cs = f1/r;
            sn = g1/r;
            for (int i=1; i<=count; ++i) {
                r *= safeMax2;
            }
        } else if (scale<=safeMin2) {
            int count = 0;
            do {
                ++count;
                f1 *= safeMax2;
                g1 *= safeMax2;
                scale = max(abs(f1), abs(g1));
            } while (scale<=safeMin2);
            r = sqrt(pow(f1,2) + pow(g1,2));
            cs = f1/r;
            sn = g1/r;
            for (int i=1; i<=count; ++i) {
                r *= safeMin2;
            }
        } else {
            r = sqrt(pow(f1,2) + pow(g1,2));
            cs = f1/r;
            sn = g1/r;
        }
        if (abs(f)>abs(g) && cs<Zero) {
            cs = -cs;
            sn = -sn;
            r = -r;
        }
    }
}

template <typename T>
T
absInf(const std::complex<T> &z)
{
    using std::abs;
    using std::conj;
    using std::imag;
    using std::max;

    return max(abs(real(z)), abs(imag(z)));
}

template <typename T>
T
absSq(const std::complex<T> &z)
{
    using cxxblas::pow;

    return pow(real(z),2) + pow(imag(z),2);
}

//
//  Complex variant
//
template <typename T, typename PT>
typename RestrictTo<IsComplex<T>::value
                 && IsReal<PT>::value,
         void>::Type
lartg_impl(const T &f, const T &g, PT &cs, T &sn, T &r)
{
    using std::conj;
    using std::real;
    using std::sqrt;

    typedef int  IndexType;

    const T   CZero(0);
    const PT  Zero(0), One(1), Two(2);

    PT   eps, safeMin, safeMin2, safeMax2;

    safeMin       = lamch<PT>(SafeMin);
    eps           = lamch<PT>(Eps);
    const int exp = explicit_cast<PT,int>(log(safeMin/eps)
                                          / log(lamch<PT>(Base))
                                          / Two);
    safeMin2      = pow(lamch<PT>(Base), exp);
    safeMax2      = One / safeMin2;

    PT scale = max(absInf(f), absInf(g));
    T  fs    = f;
    T  gs    = g;

    IndexType count = 0;
    if (scale >= safeMax2) {
        do {
            ++count;
            fs *= safeMin2;
            gs *= safeMin2;
            scale *= safeMin2;
        } while (scale >= safeMax2);
    } else if (scale<=safeMin2) {
        if (g == CZero) {
            cs = One;
            sn = CZero;
            r = f;
            return;
        }
        do {
            --count;
            fs *= safeMax2;
            gs *= safeMax2;
            scale *= safeMax2;
        } while (scale <= safeMin2);
    }
    const PT f2 = absSq(fs);
    const PT g2 = absSq(gs);

    if (f2 <= max(g2, One)*safeMin) {
//
//      This is a rare case: F is very small.
//
        if (f == CZero) {
            cs = Zero;
            r  = lapy2(real(g), imag(g));
//          Do complex/real division explicitly with two real divisions
            PT d  = lapy2(real(gs), imag(gs));
            sn    = T(real(gs) / d, -imag(gs) / d);
            return;
        }
        PT f2s = lapy2(real(fs), imag(fs));
//      G2 and G2S are accurate
//      G2 is at least SAFMIN, and G2S is at least SAFMN2
        PT g2s = sqrt(g2);
//      Error in CS from underflow in F2S is at most
//      UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS
//      If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN,
//      and so CS .lt. sqrt(SAFMIN)
//      If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN
//      and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS)
//      Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S
        cs = f2s / g2s;
//      Make sure abs(FF) = 1
//      Do complex/real division explicitly with 2 real divisions
        T ff;
        if (absInf(f) > One) {
            PT d = lapy2(real(f), imag(f));
            ff   = T(real(f)/d, imag(f)/d);
        } else {
            PT dr = safeMax2*real(f);
            PT di = safeMax2*imag(f);
            PT d  = lapy2(dr, di);
            ff    = T(dr/d, di/d);
        }
        sn = ff * T(real(gs)/g2s, -imag(gs)/g2s);
        r  = cs*f + sn*g;
    } else {
//
//      This is the most common case.
//      Neither F2 nor F2/G2 are less than SAFMIN
//      F2S cannot overflow, and it is accurate
//
        PT f2s = sqrt(One + g2/f2);
//      Do the F2S(real)*FS(complex) multiply with two real multiplies
        r  = T(f2s*real(fs), f2s*imag(fs));
        cs = One / f2s;
        PT d  = f2 + g2;
//      Do complex/real division explicitly with two real divisions
        sn =  T(real(r)/d, imag(r)/d);
        sn *= conj(gs);
        if (count != IndexType(0)) {
            if (count > IndexType(0)) {
                for (IndexType i=1; i<=count; ++i) {
                    r *= safeMax2;
                }
            } else {
                for (IndexType i=1; i<=-count; ++i) {
                    r *= safeMin2;
                }
            }
        }
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//
//  Real/complex variant
//
template <typename T, typename PT>
void
lartg_impl(const T &f, const T &g, PT &cs, T &sn, T &r)
{
    cxxlapack::lartg(f, g, cs, sn, r);
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//
//  Real variant
//
template <typename T>
typename RestrictTo<IsReal<T>::value,
         void>::Type
lartg(const T &f, const T &g, T &cs, T &sn, T &r)
{
    LAPACK_DEBUG_OUT("lartg");

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    T   _cs = cs;
    T   _sn = sn;
    T   _r  = r;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::lartg_impl(f, g, cs, sn, r);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    external::lartg_impl(f, g, _cs, _sn, _r);

    bool failed = false;
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

    if (! isIdentical(r, _r, " r", "_r")) {
        std::cerr << "CXXLAPACK:  r = " << r << std::endl;
        std::cerr << "F77LAPACK: _r = " << _r << std::endl;
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
template <typename T, typename PT>
typename RestrictTo<IsComplex<T>::value
                 && IsReal<PT>::value,
         void>::Type
lartg(const T &f, const T &g, PT &cs, T &sn, T &r)
{
    LAPACK_DEBUG_OUT("lartg (complex)");

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    PT  _cs = cs;
    T   _sn = sn;
    T   _r  = r;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::lartg_impl(f, g, cs, sn, r);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    external::lartg_impl(f, g, _cs, _sn, _r);

    bool failed = false;
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

    if (! isIdentical(r, _r, " r", "_r")) {
        std::cerr << "CXXLAPACK:  r = " << r << std::endl;
        std::cerr << "F77LAPACK: _r = " << _r << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LARTG_TCC
