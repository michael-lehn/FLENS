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
       SUBROUTINE DRSCL( N, SA, SX, INCX )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_EIG_RSCL_TCC
#define FLENS_LAPACK_EIG_RSCL_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
template <typename SA, typename VSX>
void
rscl_generic(const SA &sa, DenseVector<VSX> &sx)
{
    using std::abs;

    typedef typename DenseVector<VSX>::ElementType  ElementType;
    typedef typename DenseVector<VSX>::IndexType    IndexType;

    const IndexType n = sx.length();

    const ElementType Zero(0), One(1);
//
//  Quick return if possible
//
    if (n==0) {
        return;
    }
//
//  Get machine parameters
//
    ElementType smallNum = lamch<ElementType>(SafeMin);
    ElementType bigNum = One / smallNum;
    labad(smallNum, bigNum);
//
//  Initialize the denominator to SA and the numerator to 1.
//
    ElementType cDen = sa;
    ElementType cNum = One;

    bool done = false;
    do {
        const ElementType cDen1 = cDen*smallNum;
        const ElementType cNum1 = cNum / bigNum;

        ElementType mul;

        if (abs(cDen1)>abs(cNum) && cNum!=Zero) {
//
//          Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.
//
            mul = smallNum;
            done = false;
            cDen = cDen1;
        } else if (abs(cNum1)>abs(cDen)) {
//
//          Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.
//
            mul = bigNum;
            done = false;
            cNum = cNum1;
        } else {
//
//          Multiply X by CNUM / CDEN and return.
//
            mul = cNum / cDen;
            done = true;
        }
//
//      Scale the vector X by MUL
//
        sx *= mul;

    } while (!done);
}

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename SA, typename VSX>
void
rscl(const SA &sa, DenseVector<VSX> &sx)
{
    typedef typename DenseVector<VSX>::IndexType  IndexType;

    cxxlapack::rscl<IndexType>(sx.length(), sa, sx.data(), sx.stride());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
template <typename SA, typename VSX>
void
rscl(const SA &sa, DenseVector<VSX> &sx)
{
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(sx.firstIndex()==1);
    ASSERT(sx.inc()>0);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename DenseVector<VSX>::NoView  sx_org = sx;
#   endif

//
//  Call implementation
//
    rscl_generic(sa, sx);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename DenseVector<VSX>::NoView  sx_generic = sx;

//
//  restore output arguments
//
    sx = sx_org;

//
//  Compare generic results with results from the native implementation
//
    external::rscl(sa, sx);

    bool failed = false;
    if (! isIdentical(sx_generic, sx, "sx_generic", "sx")) {
        std::cerr << "CXXLAPACK: sx_generic = " << sx_generic << std::endl;
        std::cerr << "F77LAPACK: sx = " << sx << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename SA, typename VSX>
void
rscl(const SA &sa, VSX &&sx)
{
    CHECKPOINT_ENTER;
    rscl(sa, sx);
    CHECKPOINT_LEAVE;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_RSCL_TCC
