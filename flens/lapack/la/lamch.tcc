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
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
 *
 *  -- LAPACK auxiliary routine (version 3.3.0) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     Based on LAPACK DLAMCH but with Fortran 95 query functions
 *     See: http://www.cs.utk.edu/~luszczek/lapack/lamch.html
 *     and
 *   http://www.netlib.org/lapack-dev/lapack-coding/program-style.html#id2537289
 *     July 2010
 */

#ifndef FLENS_LAPACK_LA_LAMCH_TCC
#define FLENS_LAPACK_LA_LAMCH_TCC 1

#include <limits>

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename T>
T
lamch_impl(MachineParameter machineParameter)
{
//
//  Assume rounding, not chopping. Always.
//
    const T eps = std::numeric_limits<T>::epsilon() * T(0.5);

    if (machineParameter==Eps) {
        return eps;
    } else if (machineParameter==SafeMin) {
        T safeMin = std::numeric_limits<T>::min();
        const T small = T(1) / std::numeric_limits<T>::max();
        if (small>=safeMin) {
//
//          Use SMALL plus a bit, to avoid the possibility of rounding
//          causing overflow when computing  1/sfmin.
//
            safeMin = small*(T(1)+eps);
        }
        return safeMin;
    } else if (machineParameter==Base) {
        return std::numeric_limits<T>::radix;
    } else if (machineParameter==Precision) {
        return eps*std::numeric_limits<T>::radix;
    } else if (machineParameter==Mantissa) {
        return std::numeric_limits<T>::digits;
    } else if (machineParameter==Rounding) {
        return T(1);
    } else if (machineParameter==UnderflowExp) {
        return std::numeric_limits<T>::min_exponent;
    } else if (machineParameter==UnderflowThreshold) {
        return std::numeric_limits<T>::min();
    } else if (machineParameter==OverflowExp) {
        return std::numeric_limits<T>::max_exponent;
    } else if (machineParameter==OverflowThreshold) {
        return std::numeric_limits<T>::max();
    }
    ASSERT(0);
    return T(0);
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename T>
T
lamch_impl(MachineParameter machineParameter)
{
    return cxxlapack::lamch<T>(char(machineParameter));
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename T>
typename RestrictTo<IsNotComplex<T>::value, T>::Type
lamch(MachineParameter machineParameter)
{
//
//  Call implementation
//
    const T cxx_result = LAPACK_SELECT::lamch_impl<T>(machineParameter);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    const T f77_result = external::lamch_impl<T>(machineParameter);

    if (! isIdentical(cxx_result, f77_result, "cxx_result", "f77_result")) {
        ASSERT(0);
    }
#   endif

    return cxx_result;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAMCH_TCC
