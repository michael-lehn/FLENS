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

#ifndef FLENS_BLAS_CLOSURES_TWEAKS_R2_H
#define FLENS_BLAS_CLOSURES_TWEAKS_R2_H 1

#include <flens/blas/closures/tweaks/defaulteval.h>
#include <flens/blas/closures/tweaks/r.h>
#include <flens/blas/operators/operators.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//
//  Matrix closure of form:  a1*VX1*VY1^H + a2*VY2*VX2^H
//
template <typename VX1, typename VY1, typename VY2, typename VX2>
using MatrixClosureR2C =

    MatrixClosure<OpAdd,
                  MatrixClosureRC<VX1, VY1>,
                  MatrixClosureRC<VY2, VX2>
                 >;


//
//  Matrix closure of form:  a1*VX1*VY1^T + a2*VY2*VX2^T
//
template <typename VX1, typename VY1, typename VY2, typename VX2>
using MatrixClosureR2U =

    MatrixClosure<OpAdd,
                  MatrixClosureRU<VX1, VY1>,
                  MatrixClosureRU<VY2, VX2>
                 >;


//-- HermitianMatrix -----------------------------------------------------------
//
//  A += alpha*x*conjTrans(y) + conjugate(alpha)*y*conjTrans(x)
//
template <typename ALPHA,
          typename VX1, typename VY1,
          typename VY2, typename VX2,
          typename MA>
    typename RestrictTo<DefaultEval<MatrixClosureR2C<VX1, VY1, VY2, VX2>>::value
                     && IsVector<VX1>::value
                     && IsVector<VY1>::value
                     && IsVector<VY2>::value
                     && IsVector<VX2>::value
                     && IsHermitianMatrix<MA>::value,
             void>::Type
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureR2C<VX1, VY1, VY2, VX2> &xyh_yxh,
         MA &A);


//-- SymmetricMatrix -----------------------------------------------------------
//
//  A += alpha*x*transpose(y) + alpha*y*transpose(x)
//
template <typename ALPHA,
          typename VX1, typename VY1,
          typename VY2, typename VX2,
          typename MA>
    typename RestrictTo<DefaultEval<MatrixClosureR2U<VX1, VY1, VY2, VX2>>::value
                     && IsVector<VX1>::value
                     && IsVector<VY1>::value
                     && IsVector<VY2>::value
                     && IsVector<VX2>::value
                     && IsSymmetricMatrix<MA>::value,
             void>::Type
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureR2U<VX1, VY1, VY2, VX2> &xyt_yxt,
         MA &A);


} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_TWEAKS_R2_H
