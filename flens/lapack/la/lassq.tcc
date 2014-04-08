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
       SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LASSQ_TCC
#define FLENS_LAPACK_LA_LASSQ_TCC 1

#include <cmath>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename VX, typename T>
void
lassq_impl(const DenseVector<VX> &x, T &scale, T &sumsq)
{
    using std::abs;
    using cxxblas::pow;

    typedef typename DenseVector<VX>::IndexType IndexType;

    const IndexType n = x.length();

    if (n>0) {
        for (IndexType i=1; i<=n; ++i) {
            if (x(i)!=T(0)) {
                const T absXi = abs(x(i));
                if (scale<absXi) {
                    sumsq = 1 + sumsq*pow(scale/absXi, 2);
                    scale = absXi;
                } else {
                    sumsq += pow(absXi/scale, 2);
                }
            }
        }
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename VX, typename T>
void
lassq_impl(const DenseVector<VX> &x, T &scale, T &sumsq)
{
    cxxlapack::lassq(x.length(), x.data(), x.stride(), scale, sumsq);
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename VX, typename T>
void
lassq(const DenseVector<VX> &x, T &scale, T &sumsq)
{
    LAPACK_DEBUG_OUT("lassq");

//
//  Test the input parameters
//
    ASSERT(x.inc()>0);

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    T _scale = scale;
    T _sumsq = sumsq;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::lassq_impl(x, scale, sumsq);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    external::lassq_impl(x, _scale, _sumsq);

    bool failed = false;
    if (! isIdentical(scale, _scale, " scale", "_scale")) {
        std::cerr << "CXXLAPACK:  scale = " << scale << std::endl;
        std::cerr << "F77LAPACK: _scale = " << _scale << std::endl;
        failed = true;
    }

    if (! isIdentical(sumsq, _sumsq, " sumsq", "_sumsq")) {
        std::cerr << "CXXLAPACK:  sumsq = " << sumsq << std::endl;
        std::cerr << "F77LAPACK: _sumsq = " << _sumsq << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    } else {
//        std::cerr << "passed: lassq.tcc" << std::endl;
    }

#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LASSQ_TCC
