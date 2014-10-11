/*
 *   Copyright (c) 2012, Michael Lehn, Klaus Pototzky
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
       SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
       SUBROUTINE ZSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_SY_TRF_TCC
#define FLENS_LAPACK_SY_TRF_TCC 1

#include <cxxstd/algorithm.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (sy)trf [real variant] ----------------------------------------

template <typename MA, typename VP, typename VWORK>
typename SyMatrix<MA>::IndexType
trf_impl(SyMatrix<MA> &A, DenseVector<VP> &piv, DenseVector<VWORK> &work)
{
    typedef typename SyMatrix<MA>::IndexType   IndexType;
    typedef typename SyMatrix<MA>::ElementType ElementType;

    if (work.length()==0) {
        ElementType WORK;
        IndexType   LWORK = -1;

        cxxlapack::sytrf<IndexType>(getF77Char(A.upLo()),
                                    A.dim(),
                                    A.data(), A.leadingDimension(),
                                    piv.data(),
                                    &WORK,
                                    LWORK);
        work.resize(cxxblas::real(WORK));
    }

    return cxxlapack::sytrf<IndexType>(getF77Char(A.upLo()),
                                       A.dim(),
                                       A.data(), A.leadingDimension(),
                                       piv.data(),
                                       work.data(),
                                       work.length());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

#ifdef USE_CXXLAPACK

//-- (sy)trf [real and complex variant] --------------------------------

template <typename MA, typename VPIV, typename VWORK>
typename RestrictTo<IsSyMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsDenseVector<VWORK>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
trf(MA &&A, VPIV &&piv, VWORK &&work)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

    if (piv.length()<A.dim()) {
        piv.resize(A.dim());
    }

#   ifndef NDEBUG
//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT((piv.inc()>0 && piv.firstIndex()==1)
        || (piv.inc()<0 && piv.firstIndex()==A.numRows()));
#   endif

//
//  Call implementation
//
    IndexType info = external::trf_impl(A, piv, work);

    return info;
}


//-- (sy)trf [real and complex variant] ----------------------------------------

template <typename MA, typename VPIV>
typename RestrictTo<IsSyMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
trf(MA &&A, VPIV &&piv)
{
    typedef typename RemoveRef<MA>::Type::Vector WorkVector;

    WorkVector  work;

    return trf(A, piv, work);
}

#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_SY_TRF_TCC
