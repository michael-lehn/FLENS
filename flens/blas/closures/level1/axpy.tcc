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

#ifndef FLENS_BLAS_CLOSURES_LEVEL1_AXPY_TCC
#define FLENS_BLAS_CLOSURES_LEVEL1_AXPY_TCC 1

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
//  y += alpha*conjugate(x)
//
template <typename ALPHA, typename VL, typename VR, typename VY>
typename RestrictTo<VCDefaultEval<OpConj, VL, VR>::value
                 && IsVector<VL>::value
                 && IsVector<VR>::value,
         void>::Type
axpy(const ALPHA &alpha, const VectorClosure<OpConj, VL, VR> &x, Vector<VY> &y)
{
//
//  No need to add another log-entry as we simply pass-through to the
//  BLAS implementation
//
    acxpy(alpha, x.left(), y.impl());
}

//
//  y += alpha*(x1+x2)
//
#ifndef FLENS_DEBUG_CLOSURES

template <typename ALPHA, typename VL, typename VR, typename VY>
typename RestrictTo<VCDefaultEval<OpAdd, VL, VR>::value
                 && IsVector<VL>::value
                 && IsVector<VR>::value,
         void>::Type
axpy(const ALPHA &, const VectorClosure<OpAdd, VL, VR> &, Vector<VY> &)
{
//
//  Note that y += x1+x2  or  y -= x1+x2 is equivalent to y = y + (x1+x2) or
//  y = y - (x1+x2) respectively.  Hence, the result of x1+x2 has to be
//  computed *before* we can update y.  This means that a temporary vectors is
//  needed to hold the result of x1+x2.  We only allow this if the
//  FLENS_DEBUG_CLOSURES macro is defined such that a user has a chance to
//  optimize his/her expressions.
//
    ASSERT(0);
}

#else

template <typename ALPHA, typename VL, typename VR, typename VY>
typename RestrictTo<VCDefaultEval<OpAdd, VL, VR>::value
                 && IsVector<VL>::value
                 && IsVector<VR>::value,
         void>::Type
axpy(const ALPHA &alpha, const VectorClosure<OpAdd, VL, VR> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_AXPY(alpha, x, y);
    typedef VectorClosure<OpAdd, VL, VR>  VC;

//
//  Compute the result of closure x = (x.left()+x.right()) first and store
//  it in tmp.
//
    FLENS_BLASLOG_TMP_TRON;
    typename Result<VC>::Type tmp = x;
    FLENS_BLASLOG_TMP_TROFF;
//
//  Update y with tmp, i.e. compute y = y + alpha*tmp
//
    axpy(alpha, tmp, y.impl());

    FLENS_BLASLOG_TMP_REMOVE(tmp, x);
    FLENS_BLASLOG_END;
}

#endif

//------------------------------------------------------------------------------
//
//  y += alpha*(x1-x2)
//
#ifndef FLENS_DEBUG_CLOSURES

template <typename ALPHA, typename VL, typename VR, typename VY>
typename RestrictTo<VCDefaultEval<OpSub, VL, VR>::value
                 && IsVector<VL>::value
                 && IsVector<VR>::value,
         void>::Type
axpy(const ALPHA &, const VectorClosure<OpSub, VL, VR> &, Vector<VY> &)
{
//
//  Note that y += x1-x2  or  y -= x1-x2 is equivalent to y = y + (x1-x2) or
//  y = y - (x1-x2) respectively.  Hence, the result of x1-x2 has to be
//  computed *before* we can update y.  This means that a temporary vectors is
//  needed to hold the result of x1-x2.  We only allow this if the
//  FLENS_DEBUG_CLOSURES macro is defined such that a user has a chance to
//  optimize his/her expressions.
//
    ASSERT(0);
}

#else

template <typename ALPHA, typename VL, typename VR, typename VY>
typename RestrictTo<VCDefaultEval<OpSub, VL, VR>::value
                 && IsVector<VL>::value
                 && IsVector<VR>::value,
         void>::Type
axpy(const ALPHA &alpha, const VectorClosure<OpSub, VL, VR> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_AXPY(alpha, x, y);
    typedef VectorClosure<OpSub, VL, VR>  VC;

//
//  Compute the result of closure x = (x.left()-x.right()) first and store
//  it in tmp.
//
    FLENS_BLASLOG_TMP_TRON;
    typename Result<VC>::Type tmp = x;
    FLENS_BLASLOG_TMP_TROFF;
//
//  Update y with tmp, i.e. compute y = y + alpha*tmp
//
    axpy(alpha, tmp, y.impl());

    FLENS_BLASLOG_TMP_REMOVE(tmp, x);
    FLENS_BLASLOG_END;
}

#endif

//------------------------------------------------------------------------------
//
//  y += scalar*x
//
//  We evalute this with a scalSwitch
//  case 1: x is no closure
//  case 2: x is a closure

//
//  case 1: x is no closure
//
template <typename ALPHA, typename VX, typename VY>
void
axpySwitch(const ALPHA &alpha, const Vector<VX> &x, Vector<VY> &y)
{
//
//  No need to add another log-entry as we simply pass-through to the
//  BLAS implementation
//
    axpy(alpha, x.impl(), y.impl());
}

//
//  case 2: (a) x is a conjugate-closure
//

template <typename ALPHA, typename L, typename R, typename VY>
void
axpySwitch(const ALPHA &alpha, const VectorClosure<OpConj, L, R> &x,
           Vector<VY> &y)
{
//
//  No need to add another log-entry as we simply pass-through to the
//  BLAS implementation
//
    acxpy(alpha, x.left(), y.impl());
}

//
//  case 2: (b) x is a closure (everything else)
//
#ifndef FLENS_DEBUG_CLOSURES

template <typename ALPHA, typename Op, typename L, typename R, typename VY>
void
axpySwitch(const ALPHA &, const VectorClosure<Op, L, R> &, Vector<VY> &)
{
//  Creation of temporary not allowed

    ASSERT(0);
}

#else

template <typename ALPHA, typename Op, typename L, typename R, typename VY>
void
axpySwitch(const ALPHA &alpha, const VectorClosure<Op, L, R> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_AXPY(alpha, x, y);
    typedef VectorClosure<Op, L, R>  VC;

//
//  Compute the result of closure x and store it in tmp.
//
    FLENS_BLASLOG_TMP_TRON;
    typename Result<VC>::Type tmp = x;
    FLENS_BLASLOG_TMP_TROFF;
//
//  Update y with tmp, i.e. compute y = y + alpha*tmp
//
    axpy(alpha, tmp, y.impl());

    FLENS_BLASLOG_TMP_REMOVE(tmp, x);
    FLENS_BLASLOG_END;
}

#endif

//
//  entry point for switch
//
template <typename ALPHA, typename SV, typename VX, typename VY>
typename RestrictTo<VCDefaultEval<OpMult, SV, VX>::value
                 && IsScalarValue<SV>::value
                 && IsVector<VX>::value,
         void>::Type
axpy(const ALPHA &alpha, const VectorClosure<OpMult, SV, VX> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_AXPY(alpha, x, y);
//
//  In the entry point we either have y += x or y -= x
//
    ASSERT(alpha==ALPHA(1) || alpha==ALPHA(-1));
//
//  Switch: 1) If x is a closure we need a temporary otherwise
//          2) call BLAS implementation directly.
//
    axpySwitch(alpha*x.left().value(), x.right(), y.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  y += x/scalar
//
//  We evalute this with a rscalSwitch
//  case 1: x is no closure
//  case 2: x is a closure

//
//  case 1: x is no closure
//
template <typename ALPHA, typename VX, typename VY>
void
raxpySwitch(const ALPHA &alpha, const Vector<VX> &x, Vector<VY> &y)
{
//
//  No need to add another log-entry as we simply pass-through to the
//  BLAS implementation
//
    raxpy(alpha, x.impl(), y.impl());
}

//
//  case 2: (a) x = conjugate(x)/a
//
template <typename ALPHA, typename L, typename R, typename VY>
void
raxpySwitch(const ALPHA &alpha, const VectorClosure<OpConj, L, R> &x,
            Vector<VY> &y)
{
//  No need to add another log-entry as we simply pass-through to the
//  BLAS implementation
//
    racxpy(alpha, x.left(), y.impl());
}

//
//  case 2: (b) x is any other closure
//
#ifndef FLENS_DEBUG_CLOSURES

template <typename ALPHA, typename Op, typename L, typename R, typename VY>
void
raxpySwitch(const ALPHA &, const VectorClosure<Op, L, R> &, Vector<VY> &)
{
//
//  Creation of temporary not allowed
//
    ASSERT(0);
}

#else

template <typename ALPHA, typename Op, typename L, typename R, typename VY>
void
raxpySwitch(const ALPHA &alpha, const VectorClosure<Op, L, R> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_RAXPY(alpha, x, y);
    typedef VectorClosure<Op, L, R>  VC;

//
//  Compute the result of closure x and store it in tmp.
//
    FLENS_BLASLOG_TMP_TRON;
    typename Result<VC>::Type tmp = x;
    FLENS_BLASLOG_TMP_TROFF;
//
//  Update y with tmp, i.e. compute y = y + tmp/alpha
//
    raxpy(alpha, tmp, y.impl());

    FLENS_BLASLOG_TMP_REMOVE(tmp, x);
    FLENS_BLASLOG_END;
}

#endif

//
//  entry point for switch
//
template <typename ALPHA, typename VX, typename SV, typename VY>
typename RestrictTo<VCDefaultEval<OpDiv, VX, SV>::value
                 && IsScalarValue<SV>::value
                 && IsVector<VX>::value,
         void>::Type
axpy(const ALPHA &alpha, const VectorClosure<OpDiv, VX, SV> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_AXPY(alpha, x, y);
//
//  In the entry point we either have y += x/scalar or y -= x/scalar
//
    ASSERT(alpha==ALPHA(1) || alpha==ALPHA(-1));
//
//  Switch: 1) If x is a closure we need a temporary otherwise
//          2) call BLAS implementation directly.
//
    raxpySwitch(alpha*x.right().value(), x.left(), y.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
// y += A*x
//
template <typename ALPHA, typename ML, typename VR, typename VY>
typename RestrictTo<VCDefaultEval<OpMult, ML, VR>::value
                 && IsMatrix<ML>::value
                 && IsVector<VR>::value,
         void>::Type
axpy(const ALPHA &alpha,
     const VectorClosure<OpMult, ML, VR> &Ax, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_AXPY(alpha, Ax, y);
//
//  In the entry point we either have y += A*x or y -= A*x
//
    ASSERT(alpha==ALPHA(1) || alpha==ALPHA(-1));

//
//  If A is a closure then prune arbitrary many OpTrans/OpConj
//
    typedef typename PruneConjTrans<ML>::Remainder MA;

    Transpose trans = PruneConjTrans<ML>::trans;
    const MA  &A    = PruneConjTrans<ML>::remainder(Ax.left());
//
//  If x is a closure it gets evaluated.  In this case a temporary gets
//  created.  Otherwise we only keep a reference
//
    FLENS_BLASLOG_TMP_TRON;
    const typename Result<VR>::Type  &x = Ax.right();
    FLENS_BLASLOG_TMP_TROFF;
//
//  Call mv implementation
//
    typename VY::Impl::ElementType  One(1);
    mvSwitch(trans, alpha, A, x, One, y.impl());

//
//  If a temporary was created and registered before we now unregister it
//
#   ifdef FLENS_DEBUG_CLOSURES
    if (!IsSame<VR, typename Result<VR>::Type>::value) {
        FLENS_BLASLOG_TMP_REMOVE(x, Ax.right());
    }
#   endif

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
// y += x1 % x2
//
template <typename ALPHA, typename VL, typename VR, typename VY>
typename RestrictTo<VCDefaultEval<OpCross, VL, VR>::value
                 && IsVector<VL>::value
                 && IsVector<VR>::value,
         void>::Type
axpy(const ALPHA &alpha, const VectorClosure<OpCross, VL, VR> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_AXPY(alpha, x, y);
//
// Computes the result of the cross closure x first and stores it in vx.
//
    typedef VectorClosure<OpCross, VL, VR> VC;
    const typename Result<VC>::Type vx = x;

//
//  Update y with vx, i.e. compute y = y + alpha*vx
//
    axpy(alpha, vx.impl(), y.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
// y += x*A
//
template <typename ALPHA, typename VL, typename MR, typename VY>
typename RestrictTo<VCDefaultEval<OpMult, VL, MR>::value
                 && IsVector<VL>::value
                 && IsMatrix<MR>::value,
         void>::Type
axpy(const ALPHA &alpha, const VectorClosure<OpMult, VL, MR> &xA, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_AXPY(alpha, xA, y);
//
//  In the entry point we either have y += x*A or y -= x*A
//
    ASSERT(alpha==ALPHA(1) || alpha==ALPHA(-1));

//
//  If A is a closure then prune arbitrary many OpTrans/OpConj
//
    typedef typename PruneConjTrans<MR>::Remainder MA;

    Transpose trans = PruneConjTrans<MR>::trans;
    const MA  &A    = PruneConjTrans<MR>::remainder(xA.right());
//
//  If x is a closure it gets evaluated.  In this case a temporary gets
//  created.  Otherwise we only keep a reference
//
    FLENS_BLASLOG_TMP_TRON;
    const typename Result<VL>::Type  &x = xA.left();
    FLENS_BLASLOG_TMP_TROFF;
//
//  Call mv implementation
//
    typename VY::Impl::ElementType  One(1);
    mvSwitch(Transpose(Trans^trans), alpha, A, x, One, y.impl());

//
//  If a temporary was created and registered before we now unregister it
//
#   ifdef FLENS_DEBUG_CLOSURES
    if (!IsSame<VL, typename Result<VL>::Type>::value) {
        FLENS_BLASLOG_TMP_REMOVE(x, xA.left());
    }
#   endif

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  This gets called if everything else fails
//

//
//  y += <Unknown Closure>
//
template <typename ALPHA, typename Op, typename VL, typename VR, typename VY>
void
axpy(const ALPHA                      &FLENS_BLASLOG_VARDECL(alpha),
     const VectorClosure<Op, VL, VR>  &FLENS_BLASLOG_VARDECL(x),
     Vector<VY>                       &FLENS_BLASLOG_VARDECL(y))
{
    FLENS_BLASLOG_ERROR_AXPY(alpha, x, y);

    ASSERT(0);
}

//-- matrix closures -----------------------------------------------------------
//
//  In the following comments op(X) denotes  x, X^T or X^H
//

//
//  B += alpha*op(A1 + A2)
//
#ifndef FLENS_DEBUG_CLOSURES

template <typename ALPHA, typename ML, typename MR, typename MB>
typename RestrictTo<MCDefaultEval<OpAdd, ML, MR>::value
                 && IsMatrix<ML>::value
                 && IsMatrix<MR>::value,
         void>::Type
axpy(Transpose, const ALPHA &, const MatrixClosure<OpAdd, ML, MR> &,
     Matrix<MB> &)
{
//
//  Note that B += A1+A2  or  B -= A1+A2 is equivalent to B = B + (A1+A2) or
//  B = B - (A1+A2) respectively.  Hence, the result of A1+A2 has to be
//  computed *before* we can update y.  This means that a temporary vectors is
//  needed to hold the result of A1+A2.  We only allow this if the
//  FLENS_DEBUG_CLOSURES macro is defined such that a user has a chance to
//  optimize his/her expressions.
//
    ASSERT(0);
}

#else

template <typename ALPHA, typename ML, typename MR, typename MB>
typename RestrictTo<MCDefaultEval<OpAdd, ML, MR>::value
                 && IsMatrix<ML>::value
                 && IsMatrix<MR>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha,
     const MatrixClosure<OpAdd, ML, MR> &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);
    typedef MatrixClosure<OpAdd, ML, MR>  MC;

//
//  Compute the result of closure A = (A.left()+A.right()) first and store
//  it in tmp.
//
    FLENS_BLASLOG_TMP_TRON;
    typename Result<MC>::Type tmp = A;
    FLENS_BLASLOG_TMP_TROFF;
//
//  Update y with tmp, i.e. compute B = B + alpha*tmp
//
    axpy(trans, alpha, tmp, B.impl());

    FLENS_BLASLOG_TMP_REMOVE(tmp, A);
    FLENS_BLASLOG_END;
}

#endif

//------------------------------------------------------------------------------
//
//  B += alpha*op(A1 - A2)
//
#ifndef FLENS_DEBUG_CLOSURES

template <typename ALPHA, typename ML, typename MR, typename MB>
typename RestrictTo<MCDefaultEval<OpSub, ML, MR>::value
                 && IsMatrix<ML>::value
                 && IsMatrix<MR>::value,
         void>::Type
axpy(Transpose, const ALPHA &, const MatrixClosure<OpSub, ML, MR> &,
     Matrix<MB> &)
{
//
//  Note that B += A1+A2  or  B -= A1+A2 is equivalent to B = B + (A1+A2) or
//  B = B - (A1+A2) respectively.  Hence, the result of A1+A2 has to be
//  computed *before* we can update y.  This means that a temporary vectors is
//  needed to hold the result of A1+A2.  We only allow this if the
//  FLENS_DEBUG_CLOSURES macro is defined such that a user has a chance to
//  optimize his/her expressions.
//
    ASSERT(0);
}

#else

template <typename ALPHA, typename ML, typename MR, typename MB>
typename RestrictTo<MCDefaultEval<OpSub, ML, MR>::value
                 && IsMatrix<ML>::value
                 && IsMatrix<MR>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha,
     const MatrixClosure<OpSub, ML, MR> &A, Matrix<MB> &B)
{
    typedef MatrixClosure<OpSub, ML, MR>  MC;
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);

//
//  Compute the result of closure A = (A.left()+A.right()) first and store
//  it in tmp.
//
    FLENS_BLASLOG_TMP_TRON;
    typename Result<MC>::Type tmp = A;
    FLENS_BLASLOG_TMP_TROFF;
//
//  Update y with tmp, i.e. compute B = B + alpha*tmp
//
    axpy(trans, alpha, tmp, B.impl());

    FLENS_BLASLOG_TMP_REMOVE(tmp, A);
    FLENS_BLASLOG_END;
}

#endif

//------------------------------------------------------------------------------
//
//  B += scalar*op(A)
//

//  We evalute this with a scalSwitch
//  case 1: A is no closure
//  case 2: A is a closure

//
//  case 1: A is no closure
//
template <typename ALPHA, typename MA, typename MB>
void
axpySwitch(Transpose trans,
           const ALPHA &alpha, const Matrix<MA> &A, Matrix<MB> &B)
{
//
//  No need to add another log-entry as we simply pass-through to the
//  BLAS implementation
//
    axpy(trans, alpha, A.impl(), B.impl());
}

//
//  case 2: A is a closure
//

//
//  case 2 (a): A is a closure (= trans)
//
template <typename ALPHA, typename MA, typename MB>
void
axpySwitch(Transpose trans,
           const ALPHA &alpha, const MatrixClosureOpTrans<MA> &A, Matrix<MB> &B)
{
//
//  No need to add another log-entry as we simply pass-through to the
//  BLAS implementation
//
    trans = Transpose(trans^Trans);
    axpy(trans, alpha, A.left(), B.impl());
}

//
//  case 2 (b): A is a closure (= conj)
//
template <typename ALPHA, typename MA, typename MB>
void
axpySwitch(Transpose trans,
           const ALPHA &alpha, const MatrixClosureOpConj<MA> &A, Matrix<MB> &B)
{
//
//  No need to add another log-entry as we simply pass-through to the
//  BLAS implementation
//
    trans = Transpose(trans^Conj);
    axpy(trans, alpha, A.left(), B.impl());
}

//
//  case 2 (c): A is a closure (the rest)
//
#ifndef FLENS_DEBUG_CLOSURES

template <typename ALPHA, typename Op, typename L, typename R, typename VY>
void
axpySwitch(Transpose, const ALPHA &, const MatrixClosure<Op, L, R> &,
           Matrix<VY> &)
{
//
//  Creation of temporary not allowed
//
    ASSERT(0);
}

#else

template <typename ALPHA, typename Op, typename L, typename R, typename VY>
void
axpySwitch(Transpose trans,
           const ALPHA &alpha, const MatrixClosure<Op, L, R> &A, Matrix<VY> &B)
{
    typedef MatrixClosure<Op, L, R>  MC;

    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);

//
//  Compute the result of closure A and store it in tmp.
//
    FLENS_BLASLOG_TMP_TRON;
    typename Result<MC>::Type tmp = A;
    FLENS_BLASLOG_TMP_TROFF;
//
//  Update B with tmp, i.e. compute B = B + alpha*tmp
//
    axpy(trans, alpha, tmp, B.impl());

    FLENS_BLASLOG_TMP_REMOVE(tmp, A);
    FLENS_BLASLOG_END;
}

#   endif

//
//  entry point for switch
//
template <typename ALPHA, typename SV, typename MA, typename MB>
typename RestrictTo<MCDefaultEval<OpMult, SV, MA>::value
                 && IsScalarValue<SV>::value
                 && IsMatrix<MA>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha,
     const MatrixClosure<OpMult, SV, MA> &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);
//
//  In the entry point we either have B += alpha*A or B -= alpha*A
//
    ASSERT(alpha==ALPHA(1) || alpha==ALPHA(-1));
//
//  Switch: 1) If A is a closure we need a temporary otherwise
//          2) call BLAS implementation directly.
//

    axpySwitch(trans, alpha*A.left().value(), A.right(), B.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  B += op(A)/scalar
//

//  We evalute this with a scalSwitch
//  case 1: A is no closure
//  case 2: A is a closure

//
//  case 1: A is no closure
//
template <typename ALPHA, typename MA, typename MB>
void
raxpySwitch(Transpose trans,
            const ALPHA &alpha, const Matrix<MA> &A, Matrix<MB> &B)
{
//
//  No need to add another log-entry as we simply pass-through to the
//  BLAS implementation
//
    raxpy(trans, alpha, A.impl(), B.impl());
}

//
//  case 2: A is a closure
//
#ifndef FLENS_DEBUG_CLOSURES

template <typename ALPHA, typename Op, typename L, typename R, typename VY>
void
raxpySwitch(Transpose, const ALPHA &, const MatrixClosure<Op, L, R> &,
            Matrix<VY> &)
{
//
//  Creation of temporary not allowed
//
    ASSERT(0);
}

#else

template <typename ALPHA, typename Op, typename L, typename R, typename VY>
void
raxpySwitch(Transpose trans,
            const ALPHA &alpha, const MatrixClosure<Op, L, R> &A, Matrix<VY> &B)
{
    typedef MatrixClosure<Op, L, R>  MC;

    FLENS_BLASLOG_BEGIN_MRAXPY(trans, alpha, A, B);

//
//  Compute the result of closure A and store it in tmp.
//
    FLENS_BLASLOG_TMP_TRON;
    typename Result<MC>::Type tmp = A;
    FLENS_BLASLOG_TMP_TROFF;
//
//  Update B with tmp, i.e. compute B = B + tmp/alpha
//
    raxpy(trans, alpha, tmp, B.impl());

    FLENS_BLASLOG_TMP_REMOVE(tmp, A);
    FLENS_BLASLOG_END;
}

#endif

//
//  entry point for switch
//
template <typename ALPHA, typename MA, typename SV, typename MB>
typename RestrictTo<MCDefaultEval<OpDiv, MA, SV>::value
                 && IsMatrix<MA>::value
                 && IsScalarValue<SV>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha,
     const MatrixClosure<OpDiv, MA, SV> &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_BEGIN_MRAXPY(trans, alpha, A, B);
//
//  In the entry point we either have B += A/scalar or B -= A/scalar
//
    ASSERT(alpha==ALPHA(1) || alpha==ALPHA(-1));
//
//  Switch: 1) If A is a closure we need a temporary otherwise
//          2) call BLAS implementation directly.
//
    raxpySwitch(trans, alpha*A.right().value(), A.left(), B.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  B += op(conjugate(A))
//
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<DefaultEval<MatrixClosureOpConj<MA> >::value
                 && IsMatrix<MA>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha,
     const MatrixClosureOpConj<MA> &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);

    trans = Transpose(trans^Conj);

    axpy(trans, alpha, A.left(), B.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  B += op(A^T)
//
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<DefaultEval<MatrixClosureOpTrans<MA> >::value
                 && IsMatrix<MA>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha,
     const MatrixClosureOpTrans<MA> &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);

    trans = Transpose(trans^Trans);

    axpy(trans, alpha, A.left(), B.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  C += A*B
//
template <typename ALPHA, typename MA, typename MB, typename MC>
typename RestrictTo<MCDefaultEval<OpMult, MA, MB>::value
                 && IsMatrix<MA>::value
                 && IsMatrix<MB>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha,
     const MatrixClosure<OpMult,MA,MB> &AB, Matrix<MC> &C)
{
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, AB, C);

    typename MC::Impl::ElementType  One(1);

    if (trans==NoTrans || trans==ConjTrans) {
        mmSwitch(trans, trans, alpha, AB.left(), AB.right(), One, C.impl());
    } else {
        mmSwitch(trans, trans, alpha, AB.right(), AB.left(), One, C.impl());
    }

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  B += <Unknown Closure>
//
#ifdef FLENS_DEBUG_CLOSURES

template <typename ALPHA, typename Op, typename ML, typename MR, typename MB>
void
axpy(Transpose trans, const ALPHA &alpha,
     const MatrixClosure<Op, ML, MR> &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_ERROR_MAXPY(trans, alpha, A, B);
    ERROR_MSG("B += <Unknown Closure>");
    ASSERT(0);
}

#else

/*
template <typename ALPHA, typename Op, typename ML, typename MR, typename MB>
void
axpy(Transpose, const ALPHA &, const MatrixClosure<Op, ML, MR> &, Matrix<MB> &)
{
    ERROR_MSG("B += <Unknown Closure>");
    ASSERT(0);
}
*/

#endif


#ifndef FLENS_DEBUG_CLOSURES
//------------------------------------------------------------------------------
//
// B += Some Matrix
//

template <typename ALPHA, typename MA, typename MB>
void
axpy(Transpose          FLENS_BLASLOG_VARDECL(trans),
     const ALPHA        &FLENS_BLASLOG_VARDECL(alpha),
     const Matrix<MA>   &FLENS_BLASLOG_VARDECL(A),
     Matrix<MB>         &FLENS_BLASLOG_VARDECL(B))
{
    FLENS_BLASLOG_ERROR_MAXPY(trans, alpha, A, B);
    ERROR_MSG("B += Some Matrix");

    ASSERT(0);
}

#else

template <typename ALPHA, typename MA, typename MB>
void
axpy(Transpose trans, const ALPHA &alpha,
     const Matrix<MA> &A, Matrix<MB> &B)
{
    typedef typename MA::Impl::ElementType        TA;
    typedef GeMatrix<FullStorage<TA, ColMajor> >  RMA;

    FLENS_BLASLOG_TMP_TRON;
    RMA A_ = A.impl();
    FLENS_BLASLOG_TMP_TROFF;

    axpy(trans, alpha, A_, B.impl());

    FLENS_BLASLOG_TMP_REMOVE(A_, A);
}

#endif


} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_LEVEL1_AXPY_TCC
