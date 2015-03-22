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

#ifndef FLENS_BLAS_CLOSURES_LEVEL1_COPYSUM_TCC
#define FLENS_BLAS_CLOSURES_LEVEL1_COPYSUM_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/closures/closures.h>
#include <flens/blas/level1/level1.h>
#include <flens/blas/level2/level2.h>
#include <flens/blas/level3/level3.h>
#include <flens/typedefs.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//
//== vector closures ===========================================================
//

//
// Auxiliary function for
//     y = x1 + x2  (alpha= 1)
// or  y = x1 - x2  (alpha=-1)
//
template <typename VX1, typename ALPHA, typename VX2, typename VY>
void
copySum(const VX1 &x1, const ALPHA &alpha, const VX2 &x2, VY &y)
{
    ASSERT(alpha==ALPHA(1) || alpha==ALPHA(-1));
//
//  In debug-closure-mode we check if x2 has to stored in a temporary.
//  Otherwise an assertion gets triggered if x2 and y are identical of if x2
//  is a closure that contains y.
//
#   ifdef FLENS_DEBUG_CLOSURES
    bool tmpNeeded = DebugClosure::search(x2, y);
    typename Result<VX2>::NoView tmp;

    if (tmpNeeded) {
        tmp = x2;
        FLENS_BLASLOG_TMP_ADD(tmp);
    }
#   else
    ASSERT(!DebugClosure::search(x2, y));
#   endif

//
//  y = x1
//
    if (! DEBUGCLOSURE::identical(x1, y)) {
        blas::copy(x1, y);
    }
//
//  y += x2  or  y -= x2
//
#   ifdef FLENS_DEBUG_CLOSURES
    if (tmpNeeded) {
        blas::axpy(alpha, tmp, y);
        FLENS_BLASLOG_TMP_REMOVE(tmp, x2);
    } else {
        blas::axpy(alpha, x2, y);
    }
#   else
    blas::axpy(alpha, x2, y);
#   endif
}

//
// Specialization for
// y = beta*y + alpha*x
// y = beta*y - alpha*x
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
        const VectorClosure<OpMult, VXL2, VXR2> &x2, VY &y)
{
    ASSERT(alpha==ALPHA(1) || alpha==ALPHA(-1));

    ASSERT(!DebugClosure::search(x2, y));

    if (DEBUGCLOSURE::identical(x1.right().impl(), y)) {
         blas::axpby(alpha*x2.left().value(), x2.right(), x1.left().value(), y);
    } else {
         blas::copy(x1.right(), y);
         blas::scal(x1.left().value(), y);
         blas::axpy(alpha*x2.left().value(), x2.right(), y);
    }
}

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
        const VectorClosure<OpMult, VXL2, VectorClosureOpConj<VXR2> > &x2,
        VY &y)
{
    ASSERT(alpha==ALPHA(1) || alpha==ALPHA(-1));

    ASSERT(!DebugClosure::search(x2, y));

    if (DEBUGCLOSURE::identical(x1.right().impl(), y)) {
         blas::acxpby(alpha*x2.left().value(),
                      x2.right().left(),
                      x1.left().value(),
                      y);
    } else {
         blas::copy(x1.right(), y);
         blas::scal(x1.left().value(), y);
         blas::acxpy(alpha*x2.left().value(), x2.right().left(), y);
    }
}

//
//== matrix closures ===========================================================
//

//
// Auxiliary function for
//     B = op(A1 + A2)  (alpha= 1)
// or  B = op(A1 - A2)  (alpha=-1)
//
template <typename MA1, typename ALPHA, typename MA2, typename MB>
void
copySum(Transpose trans,
    const MA1 &A1, const ALPHA &alpha, const MA2 &A2, MB &B)
{
    ASSERT(alpha==ALPHA(1) || alpha==ALPHA(-1));
//
//  In debug-closure-mode we check if A2 has to stored in a temporary.
//  Otherwise an assertion gets triggered if A2 and B are identical of if A2
//  is a closure that contains B.
//
#   ifdef FLENS_DEBUG_CLOSURES
    bool tmpNeeded = DebugClosure::search(A2, B);
    typename Result<MA2>::NoView tmp;

    if (tmpNeeded) {
        tmp = A2;
        FLENS_BLASLOG_TMP_ADD(tmp);
    }
#   else
    ASSERT(!DebugClosure::search(A2, B));
#   endif

//
//  B = A1
//
    if (! DEBUGCLOSURE::identical(A1, B)) {
        blas::copy(trans, A1, B);
    }
//
//  B += A2  or  B -= A2
//
#   ifdef FLENS_DEBUG_CLOSURES
    if (tmpNeeded) {
        blas::axpy(trans, alpha, tmp, B);
        FLENS_BLASLOG_TMP_REMOVE(tmp, A2);
    } else {
        blas::axpy(trans, alpha, A2, B);
    }
#   else
    blas::axpy(trans, alpha, A2, B);
#   endif

}

//
// Auxiliary function for
//     B = beta1*op(A1) + beta2*op(A2)  (alpha= 1)
// or  B = beta1*op(A1) - beta2*op(A2)  (alpha= -1)
//
template <typename MAL1, typename MAR1, typename ALPHA,
          typename MAL2, typename MAR2, typename MB>
typename RestrictTo<IsScalarValue<MAL1>::value &&
                    IsMatrix<MAR1>::value &&
                    IsScalarValue<MAL2>::value &&
                    IsMatrix<MB>::value,
                    void>::Type
copySum(Transpose trans,
        const MatrixClosure<OpMult, MAL1, MAR1> &A1,
        const ALPHA &DEBUG_VAR(alpha),
        const MatrixClosure<OpMult, MAL2, MAR2> &A2, MB &B)
{

    typedef typename PruneConjTrans<MAR2>::Remainder  MAR2_;

    ASSERT(alpha==ALPHA(1) || alpha==ALPHA(-1));

    ASSERT(!DebugClosure::search(A2, B));

    Transpose trans_ = Transpose(trans^PruneConjTrans<MAR2>::trans);
    const MAR2_ &A2_  = PruneConjTrans<MAR2>::remainder(A2.right());
//
//  B = A1
//
    if (DEBUGCLOSURE::identical(A1.right().impl(), B)) {
        blas::axpby(trans_, A2.left().value(), A2_, A1.left().value(), B);
    } else {
        blas::copy(trans, A1.right(), B);
        blas::scal(A1.left().value(), B);
        blas::axpy(trans_, A2.left().value(), A2_, B);
    }

}


} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_LEVEL1_COPY_TCC
