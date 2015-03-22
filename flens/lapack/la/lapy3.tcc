/*
 *   Copyright (c) 2013, Michael Lehn
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
      DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LAPY3_TCC
#define FLENS_LAPACK_LA_LAPY3_TCC 1

#include <cxxstd/cmath.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename T>
T
lapy3_impl(const T &x, const T &y, const T &z)
{
    using std::abs;
    using flens::max;
    using std::min;
    using std::sqrt;

    const T Zero(0);

    const T xAbs = abs(x);
    const T yAbs = abs(y);
    const T zAbs = abs(z);

    const T w = max(xAbs, yAbs, zAbs);

    if (w==Zero) {
//      W can be zero for max(0,nan,0)
//      adding all three entries together will make sure
//      NaN will not disappear.
        return xAbs + yAbs + zAbs;
    } else {
        return w*sqrt(pow(xAbs/w,2) + pow(yAbs/w,2) + pow(zAbs/w,2));
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename T>
T
lapy3_impl(const T &x, const T &y, const T &z)
{
    return cxxlapack::lapy3(x, y, z);
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename T>
T
lapy3(const T &x, const T &y, const T &z)
{
    LAPACK_DEBUG_OUT("lapy3");

//
//  Call implementation
//
    const T result = LAPACK_SELECT::lapy3_impl(x, y, z);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    const T result_ = external::lapy3_impl(x, y, z);

    bool failed = false;
    if (! isIdentical(result, result_, " result", "result_")) {
        std::cerr << "CXXLAPACK:  result = " << result << std::endl;
        std::cerr << "F77LAPACK: result_ = " << result_ << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "x = " << x << std::endl;
        std::cerr << "y = " << y << std::endl;
        std::cerr << "z = " << z << std::endl;

        std::cerr << "hex(x) = " << hex(x) << std::endl;
        std::cerr << "hex(y) = " << hex(y) << std::endl;
        std::cerr << "hex(z) = " << hex(z) << std::endl;

        std::cerr << "hex(result)  = " << hex(result) << std::endl;
        std::cerr << "hex(result_) = " << hex(result_) << std::endl;
        ASSERT(0);
    }
#   endif

    return result;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAPY§TCC_
