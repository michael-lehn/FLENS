/*
 *   Copyright (c) 2007, Michael Lehn
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

#ifndef FLENS_BLAS_LEVEL3_MM_H
#define FLENS_BLAS_LEVEL3_MM_H

#include <cxxblas/cxxblas.h>
#include <flens/auxiliary/auxiliary.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>

namespace flens { namespace blas {

//== GeneralMatrix - GeneralMatrix products ====================================

    
    //-- diagmm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    typename RestrictTo<IsDiagMatrix<MA>::value
                     && IsGeMatrix<MB>::value
                     && IsGeMatrix<MC>::value,
             void>::Type
    mm(Transpose        transA,
       Transpose        transB,
       const ALPHA      &alpha,
       const MA         &A,
       const MB         &B,
       const BETA       &beta,
       MC               &&C);
    
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    typename RestrictTo<IsGeMatrix<MA>::value
                     && IsDiagMatrix<MB>::value
                     && IsGeMatrix<MC>::value,
             void>::Type
    mm(Transpose        transA,
       Transpose        transB,
       const ALPHA      &alpha,
       const MA         &A,
       const MB         &B,
       const BETA       &beta,
       MC               &&C);
    
//-- gbmm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    typename RestrictTo<IsGbMatrix<MA>::value
                     && IsGeMatrix<MB>::value
                     && IsGeMatrix<MC>::value,
             void>::Type
    mm(Transpose        transA,
       Transpose        transB,
       const ALPHA      &alpha,
       const MA         &A,
       const MB         &B,
       const BETA       &beta,
       MC               &&C);

template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    typename RestrictTo<IsGeMatrix<MA>::value
                     && IsGbMatrix<MB>::value
                     && IsGeMatrix<MC>::value,
             void>::Type
    mm(Transpose        transA,
       Transpose        transB,
       const ALPHA      &alpha,
       const MA         &A,
       const MB         &B,
       const BETA       &beta,
       MC               &&C);

//-- gemm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    typename RestrictTo<IsGeMatrix<MA>::value
                     && IsGeMatrix<MB>::value
                     && IsGeMatrix<MC>::value,
             void>::Type
    mm(Transpose        transA,
       Transpose        transB,
       const ALPHA      &alpha,
       const MA         &A,
       const MB         &B,
       const BETA       &beta,
       MC               &&C);

//== HermitianMatrix - GeneralMatrix products ==================================

//-- hbmm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    typename RestrictTo<IsHbMatrix<MA>::value
                     && IsGeMatrix<MB>::value
                     && IsGeMatrix<MC>::value,
             void>::Type
    mm(Side             side,
       const ALPHA      &alpha,
       const MA         &A,
       const MB         &B,
       const BETA       &beta,
       MC               &&C);

//-- hemm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    typename RestrictTo<IsHeMatrix<MA>::value
                     && IsGeMatrix<MB>::value
                     && IsGeMatrix<MC>::value,
             void>::Type
    mm(Side             side,
       const ALPHA      &alpha,
       const MA         &A,
       const MB         &B,
       const BETA       &beta,
       MC               &&C);

//== SymmetricMatrix - GeneralMatrix products ==================================

//-- sbmm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    typename RestrictTo<IsSbMatrix<MA>::value
                     && IsGeMatrix<MB>::value
                     && IsGeMatrix<MC>::value,
                 void>::Type
    mm(Side             side,
       const ALPHA      &alpha,
       const MA         &A,
       const MB         &B,
       const BETA       &beta,
       MC               &&C);

//-- symm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    typename RestrictTo<IsSyMatrix<MA>::value
                     && IsGeMatrix<MB>::value
                     && IsGeMatrix<MC>::value,
             void>::Type
    mm(Side             side,
       const ALPHA      &alpha,
       const MA         &A,
       const MB         &B,
       const BETA       &beta,
       MC               &&C);

//== TriangularMatrix - GeneralMatrix products =================================

//-- tbmm
template <typename ALPHA, typename MA, typename MB>
    typename RestrictTo<IsTbMatrix<MA>::value
                     && IsGeMatrix<MB>::value,
             void>::Type
    mm(Side             side,
       Transpose        transA,
       const ALPHA      &alpha,
       const MA         &A,
       MB               &&B);

//-- trmm
template <typename ALPHA, typename MA, typename MB>
    typename RestrictTo<IsTrMatrix<MA>::value
                     && IsGeMatrix<MB>::value,
             void>::Type
    mm(Side             side,
       Transpose        transA,
       const ALPHA      &alpha,
       const MA         &A,
       MB               &&B);

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL3_MM_H
