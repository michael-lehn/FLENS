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
       SUBROUTINE DGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,
      $                   INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_EIG_BAK_TCC
#define FLENS_LAPACK_EIG_BAK_TCC 1

#include <flens/aux/aux.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename IndexType, typename VSCALE, typename MV>
void
bak_generic(BALANCE::Balance            job,
            Side                        side,
            IndexType                   iLo,
            IndexType                   iHi,
            const DenseVector<VSCALE>   &scale,
            GeMatrix<MV>                &V)
{
    using namespace BALANCE;

    typedef typename GeMatrix<MV>::ElementType  T;
    const T One(1);

    const Underscore<IndexType> _;

    const IndexType n = V.numRows();
    const IndexType m = V.numCols();
//
//  Quick return if possible
//
    if (n==0) {
        return;
    }
    if (m==0) {
        return;
    }
    if (job==None) {
        return;
    }

    if (iLo!=iHi) {
//
//      Backward balance
//
        if (job==ScaleOnly || job==Both) {

            if (side==Right) {
                for (IndexType i=iLo; i<=iHi; ++i) {
                    V(i,_) *= scale(i);
                }
            }

            if (side==Left) {
                for (IndexType i=iLo; i<=iHi; ++i) {
                    V(i,_) *= One / scale(i);
                }
            }

        }
    }
//
//  Backward permutation
//
//  For  I = ILO-1 step -1 until 1,
//           IHI+1 step 1 until N do --
//
    if (job==PermuteOnly || job==Both) {
        if (side==Right) {
            for (IndexType ii=1; ii<=n; ++ii) {
                IndexType i = ii;
                if (i>=iLo && i<=iHi) {
                    continue;
                }
                if (i<iLo) {
                    i = iLo - ii;
                }
                const IndexType k = explicit_cast<T,IndexType>(scale(i));
                if (k==i) {
                    continue;
                }
                blas::swap(V(i,_), V(k,_));
            }
        }

        if (side==Left) {
            for (IndexType ii=1; ii<=n; ++ii) {
                IndexType i = ii;
                if (i>=iLo && i<=iHi) {
                    continue;
                }
                if (i<iLo) {
                    i = iLo - ii;
                }
                const IndexType k = explicit_cast<T,IndexType>(scale(i));
                if (k==i) {
                    continue;
                }
                blas::swap(V(i,_), V(k,_));
            }
        }
    }
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename IndexType, typename VSCALE, typename MV>
void
bak_native(BALANCE::Balance             job,
           Side                         side,
           IndexType                    iLo,
           IndexType                    iHi,
           const DenseVector<VSCALE>    &scale,
           GeMatrix<MV>                 &V)
{
    typedef typename GeMatrix<MV>::ElementType  T;

    const char       JOB    = getF77LapackChar<BALANCE::Balance>(job);
    const char       SIDE   = getF77LapackChar<Side>(side);
    const INTEGER    N      = V.numRows();
    const INTEGER    ILO    = iLo;
    const INTEGER    IHI    = iHi;
    const INTEGER    M      = V.numCols();
    const INTEGER    LDV    = V.leadingDimension();
    INTEGER          INFO;

    if (IsSame<T, DOUBLE>::value) {
        LAPACK_IMPL(dgebak)(&JOB,
                            &SIDE,
                            &N,
                            &ILO,
                            &IHI,
                            scale.data(),
                            &M,
                            V.data(),
                            &LDV,
                            &INFO);
    } else {
        ASSERT(0);
    }
    ASSERT(INFO==0);
}

#endif

//== public interface ==========================================================

template <typename IndexType, typename VSCALE, typename MV>
void
bak(BALANCE::Balance            job,
    Side                        side,
    IndexType                   iLo,
    IndexType                   iHi,
    const DenseVector<VSCALE>   &scale,
    GeMatrix<MV>                &V)
{
    LAPACK_DEBUG_OUT("bak");

#   ifndef NDEBUG
//
//  Test the input parameters
//
    const IndexType n = V.numRows();

    if (n>0) {
        ASSERT(1<=iLo);
        ASSERT(iLo<=iHi);
        ASSERT(iHi<=n);
    } else {
        ASSERT(iLo==1);
        ASSERT(iHi==0);
    }
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MV>::NoView   _V = V;
#   endif

//
//  Call implementation
//
    bak_generic(job, side, iLo, iHi, scale, V);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    bak_native(job, side, iLo, iHi, scale, _V);

    if (! isIdentical(V, _V, " V", "_V")) {
        std::cerr << "CXXLAPACK:  V = " << V << std::endl;
        std::cerr << "F77LAPACK: _V = " << _V << std::endl;
        ASSERT(0);
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename SCALE, typename MV>
void
bak(BALANCE::Balance    job,
    Side                side,
    IndexType           iLo,
    IndexType           iHi,
    const SCALE         &scale,
    MV                  &&V)
{
    CHECKPOINT_ENTER;
    bak(job, side, iLo, iHi, scale, V);
    CHECKPOINT_LEAVE;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_BAK_TCC
