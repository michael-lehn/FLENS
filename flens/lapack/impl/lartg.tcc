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
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_LARTG_TCC
#define FLENS_LAPACK_IMPL_LARTG_TCC 1

#include <flens/aux/aux.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename T>
void
lartg_impl(const T &f, const T &g, T &cs, T &sn, T &r)
{
    using std::abs;
    using std::max;
    using flens::pow;

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

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename T>
void
lartg_impl(const T &f, const T &g, T &cs, T &sn, T &r)
{
    cxxlapack::lartg(f, g, cs, sn, r);
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename T>
void
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

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_LARTG_TCC
