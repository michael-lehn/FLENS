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
      SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
     $                   WORK, LWORK, INFO )
      SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
    $                    WORK, LWORK, RWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_SVD_H
#define FLENS_LAPACK_IMPL_SVD_H 1

#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

namespace SVD {

    enum Job {
        All          = 'A', // = 'A': All columns are returned;
        Save         = 'S', // = 'S': Only min(m,n) columns are returned;
        Overwrite    = 'O', // = 'O': min(m,n) columns are saved in A
        None         = 'N'  // = 'N': No rows are computed.
    };

}


//== (ge)svd ===================================================================
#ifdef USE_CXXLAPACK
//
//  Real variant
//
template <typename MA, typename VS, typename MU, typename MVT, typename VWORK>
    typename RestrictTo<IsRealGeMatrix<MA>::value
                     && IsRealDenseVector<VS>::value
                     && IsRealGeMatrix<MU>::value
                     && IsRealGeMatrix<MVT>::value
                     && IsRealDenseVector<VWORK>::value,
             void>::Type
    svd(SVD::Job    jobU,
        SVD::Job    jobVT,
        MA          &&A,
        VS          &&s,
        MU          &&U,
        MVT         &&VT,
        VWORK       &&work);

//
//  Complex variant
//
template <typename MA, typename VS, typename MU, typename MVT, typename VWORK,
          typename VRWORK>
    typename RestrictTo<IsComplexGeMatrix<MA>::value
                     && IsRealDenseVector<VS>::value
                     && IsComplexGeMatrix<MU>::value
                     && IsComplexGeMatrix<MVT>::value
                     && IsComplexDenseVector<VWORK>::value
                     && IsRealDenseVector<VRWORK>::value,
             void>::Type
    svd(SVD::Job    jobU,
        SVD::Job    jobVT,
        MA          &&A,
        VS          &&s,
        MU          &&U,
        MVT         &&VT,
        VWORK       &&work,
        VRWORK      &&rwork);

//
//  Real variant with temporary workspace
//
template <typename MA, typename VS, typename MU, typename MVT>
    typename RestrictTo<IsRealGeMatrix<MA>::value
                     && IsRealDenseVector<VS>::value
                     && IsRealGeMatrix<MU>::value
                     && IsRealGeMatrix<MVT>::value,
             void>::Type
    svd(SVD::Job    jobU,
        SVD::Job    jobVT,
        MA          &&A,
        VS          &&s,
        MU          &&U,
        MVT         &&VT);



//
//  Complex variant with temporary workspace
//
template <typename MA, typename VS, typename MU, typename MVT>
    typename RestrictTo<IsComplexGeMatrix<MA>::value
                     && IsRealDenseVector<VS>::value
                     && IsComplexGeMatrix<MU>::value
                     && IsComplexGeMatrix<MVT>::value,
             void>::Type
    svd(SVD::Job    jobU,
        SVD::Job    jobVT,
        MA          &&A,
        VS          &&s,
        MU          &&U,
        MVT         &&VT);


//== workspace query ===========================================================

template <typename MA, typename VS, typename MU, typename MVT, typename VWORK>
    typename RestrictTo<IsGeMatrix<MA>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    svd_wsq(SVD::Job    jobU,
            SVD::Job    jobVT,
            MA          &&A);

#endif // USE_CXXLAPACK


} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_SVD_H
