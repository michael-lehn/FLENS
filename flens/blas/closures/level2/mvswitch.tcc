/*
 *   Copyright (c) 2012, Michael Lehn
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

#ifndef FLENS_BLAS_CLOSURES_LEVEL2_MVSWITCH_TCC
#define FLENS_BLAS_CLOSURES_LEVEL2_MVSWITCH_TCC 1

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
//  This switch evaluates closures of the form
//
//      y = beta*y + A*x
//
//  If x is a closure then it gets evaluated and a temporary gets created to
//  store the result.  For matrix A we distinguish between three cases:
//  case 1: A is no closure
//  case 2: A is a scaling closure (i.e. scale*A)
//  case 3: A is some other closure

//
//  Entry point for mvSwitch
//
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsSame<MA, typename MA::Impl>::value &&
                    IsSame<VX, typename VX::Impl>::value &&
                    IsSame<VY, typename VY::Impl>::value,
         void>::Type
mvSwitch(Transpose trans, const ALPHA &alpha, const MA &A, const VX &x,
         const BETA &beta, VY &y)
{
    ASSERT(alpha==ALPHA(1) || alpha==ALPHA(-1));
//
//  If A is a closure then prune arbitrary many OpTrans/OpConj
//
    typedef typename PruneConjTrans<MA>::Remainder RMA;

    trans = Transpose(trans^PruneConjTrans<MA>::trans);
    const RMA  &_A = PruneConjTrans<MA>::remainder(A);
//
//  If x is a closure it gets evaluated.  In this case a temporary gets
//  created.  Otherwise we only keep a reference
//
    FLENS_BLASLOG_TMP_TRON;
    const typename Result<VX>::Type  &_x = x;
    FLENS_BLASLOG_TMP_TROFF;
//
//  Call mv implementation
//
    mvCase(trans, alpha, _A, _x, beta, y);
//
//  If a temporary was created and registered before we now unregister it
//
#   ifdef FLENS_DEBUG_CLOSURES
    if (!IsSame<VX, typename Result<VX>::Type>::value) {
        FLENS_BLASLOG_TMP_REMOVE(_x, x);
    }
#   else
    const bool check = IsSame<VX, typename Result<VX>::Type>::value;
    if (!check) {
        std::cerr << "ERROR: Temporary required." << std::endl;
    }
    ASSERT(check);
#   endif
}

//
//  case 1: A is no closure
//
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<!IsClosure<MA>::value,
         void>::Type
mvCase(Transpose trans, const ALPHA &alpha, const MA &A, const VX &x,
       const BETA &beta, VY &y)
{
    mv(trans, alpha, A, x, beta, y);
}

//
//  case 2: A is closure of type scale*A
//
template <typename ALPHA, typename T, typename MA, typename VX, typename BETA,
          typename VY>
void
mvCase(Transpose trans, const ALPHA &alpha,
       const MatrixClosure<OpMult, ScalarValue<T>, MA> &scale_A,
       const VX &x, const BETA &beta, VY &y)
{
//
//  If A is a closure then prune arbitrary many OpTrans/OpConj
//
    typedef typename PruneConjTrans<MA>::Remainder  _MA;
    typedef typename Result<_MA>::Type              RMA;

    Transpose _trans = Transpose(trans^PruneConjTrans<MA>::trans);
//
//  If the remaining A is a closure it gets evaluated.  In this case
//  a temporary gets created.  Otherwise we only keep a reference
//
    FLENS_BLASLOG_TMP_TRON;
    const _MA &_A  = PruneConjTrans<MA>::remainder(scale_A.right());
    const RMA &A   = _A;
    FLENS_BLASLOG_TMP_TROFF;

    mv(_trans, alpha*scale_A.left().value(), A, x, beta, y);

//
//  If a temporary was created and registered before we now unregister it
//
#   ifdef FLENS_DEBUG_CLOSURES
    if (!IsSame<_MA, RMA>::value) {
        FLENS_BLASLOG_TMP_REMOVE(A, _A);
    }
#   else
    const bool check = IsSame<_MA, RMA>::value;
    if (!check) {
        std::cerr << "ERROR: Temporary required." << std::endl;
    }
    ASSERT(check);
#   endif
}

//
//  case 3: A is some other closure
//
template <typename ALPHA, typename Op, typename L, typename R, typename VX,
          typename BETA, typename VY>
void
mvCase(Transpose trans, const ALPHA &alpha, const MatrixClosure<Op, L, R> &A,
         const VX &x, const BETA &beta, VY &y)
{
    typedef MatrixClosure<Op, L, R>  MC;

//
//  Create (most certainly) temporary for the result of A
//
    FLENS_BLASLOG_TMP_TRON;
    typedef typename Result<MC>::Type  MA;
    const MA &_A = A;
    FLENS_BLASLOG_TMP_TROFF;

    mv(trans, alpha, _A, x, beta, y);

//
//  If a temporary was created and registered before we now unregister it
//
#   ifdef FLENS_DEBUG_CLOSURES
    if (!IsSame<MC, MA>::value) {
        FLENS_BLASLOG_TMP_REMOVE(_A, A);
    }
#   else
    const bool check = IsSame<MC, MA>::value;
    if (!check) {
        std::cerr << "ERROR: Temporary required." << std::endl;
    }
    ASSERT(check);
#   endif
}

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_LEVEL2_MVSWITCH_TCC
