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
      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_AUX_LAPY2_TCC
#define FLENS_LAPACK_AUX_LAPY2_TCC 1

#include <cmath>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename T>
T
lapy2_generic(const T &x, const T &y)
{
    using std::abs;
    using std::max;
    using std::min;
    using std::pow;
    using std::sqrt;

    const T xAbs = abs(x);
    const T yAbs = abs(y);

    const T w = max(xAbs, yAbs);
    const T z = min(xAbs, yAbs);

    if (z==T(0)) {
        return w;
    }

    return w*sqrt(T(1) + pow(z/w,2));
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename T>
T
lapy2_native(const T &x, const T &y)
{
    if (IsSame<T, DOUBLE>::value) {
        return LAPACK_IMPL(dlapy2)(&x,
                                   &y);
    } else {
        ASSERT(0);
    }
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename T>
T
lapy2(const T &x, const T &y)
{
    LAPACK_DEBUG_OUT("lapy2");

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    T _x    = x;
    T _y    = y;
#   endif

//
//  Call implementation
//
    const T result = lapy2_generic(x, y);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    const T _result = lapy2_native(_x, _y);

    bool failed = false;
    if (! isIdentical(x, _x, " x", "_x")) {
        std::cerr << "CXXLAPACK:  x = " << x << std::endl;
        std::cerr << "F77LAPACK: _x = " << _x << std::endl;
        failed = true;
    }
    if (! isIdentical(y, _y, " y", "_y")) {
        std::cerr << "CXXLAPACK:  y = " << y << std::endl;
        std::cerr << "F77LAPACK: _y = " << _y << std::endl;
        failed = true;
    }
    if (! isIdentical(result, _result, " result", "_result")) {
        std::cerr << "CXXLAPACK:  result = " << result << std::endl;
        std::cerr << "F77LAPACK: _result = " << _result << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif

    return result;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_LAPY2_TCC
