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

#ifndef FLENS_BLAS_CLOSURES_LEVEL1_COPYSUM_H
#define FLENS_BLAS_CLOSURES_LEVEL1_COPYSUM_H 1

#include <cxxblas/cxxblas.h>
#include <flens/blas/closures/tweaks/defaulteval.h>
#include <flens/blas/operators/operators.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//
//== for vector closures =======================================================
//

//
// Auxiliary function for
//     y = x1 + x2  (alpha= 1)
// or  y = x1 - x2  (alpha=-1)
//

template <typename VX1, typename ALPHA, typename VX2, typename VY>
    void
    copySum(const VX1 &x1, const ALPHA &alpha, const VX2 &x2, VY &y);

//
// Auxiliary function for
//     y = beta1*x1 + beta2*x2  (alpha= 1)
// or  y = beta1*x1 - beta2*x2  (alpha=-1)
//
template <typename VXL1, typename VXR1, typename ALPHA, 
          typename VXL2, typename VXR2, typename VY>
    typename RestrictTo<IsScalarValue<VXL1>::value &&
                        IsVector<VXR1>::value &&
                        IsScalarValue<VXL2>::value &&
                        IsVector<VXR2>::value &&
                        IsVector<VY>::value,
                        void>::Type
    copySum(const VectorClosure<OpMult, VXL1, VXR1> &x1, const ALPHA &alpha,
            const VectorClosure<OpMult, VXL2, VXR2> &x2, VY &y);
    
//
// Auxiliary function for
//     y = beta1*x1 + beta2*conj(x2) (alpha= 1)
// or  y = beta1*x1 - beta2*conj(x2) (alpha=-1)
//
template <typename VXL1, typename VXR1, typename ALPHA,
          typename VXL2, typename VXR2, typename VY>
          typename RestrictTo<IsScalarValue<VXL1>::value &&
                              IsVector<VXR1>::value &&
                              IsScalarValue<VXL2>::value &&
                              IsVector<VXR2>::value &&
                              DefaultEval<VectorClosureOpConj<VXR2> >::value &&
                              IsVector<VY>::value,
                              void>::Type
    copySum(const VectorClosure<OpMult, VXL1, VXR1> &x1, const ALPHA &alpha,
            const VectorClosure<OpMult, VXL2, VectorClosureOpConj<VXR2> > &x2, VY &y);


//
//== for matrix closures =======================================================
//

//
// Auxiliary function for
//     B = op(A1 + A2)  (alpha= 1)
// or  B = op(A1 - A2)  (alpha=-1)
//

template <typename MA1, typename ALPHA, typename MA2, typename MB>
    void
    copySum(Transpose trans,
            const MA1 &A1, const ALPHA &alpha, const MA2 &A2, MB &B);
    
//
// Auxiliary function for
//     B = beta1*A1 + beta2*A2  (alpha= 1)
// or  B = beta1*A1 - beta2*A2  (alpha= -1)
//
template <typename MAL1, typename MAR1, typename ALPHA,
          typename MAL2, typename MAR2, typename MB>
          typename RestrictTo<IsScalarValue<MAL1>::value &&
                              IsMatrix<MAR1>::value &&
                              IsScalarValue<MAL2>::value &&
                              IsMatrix<MB>::value,
                              void>::Type
    copySum(Transpose trans,
	        const MatrixClosure<OpMult, MAL1, MAR1> &A1, const ALPHA &alpha,
            const MatrixClosure<OpMult, MAL2, MAR2> &A2, MB &B);

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_LEVEL1_COPY_H
