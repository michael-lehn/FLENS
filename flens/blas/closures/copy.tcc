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

#ifndef FLENS_BLAS_CLOSURES_COPY_TCC
#define FLENS_BLAS_CLOSURES_COPY_TCC 1

#include <flens/aux/aux.h>
#include <flens/blas/closures/debugclosure.h>
#include <flens/blas/closures/mmswitch.h>
#include <flens/blas/closures/mvswitch.h>
#include <flens/blas/closures/result.h>
#include <flens/blas/closures/residual.h>
#include <flens/blas/level1/level1.h>
#include <flens/blas/level2/level2.h>
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
    blas::copy(x1, y);
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

//------------------------------------------------------------------------------
//
// y = x1 + x2
//
template <typename VL, typename VR, typename VY>
typename RestrictTo<!ClosureType<OpAdd, VL, VR>::isMatrixVectorProduct,
         void>::Type
copy(const VectorClosure<OpAdd, VL, VR> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(x, y);

    typedef typename VY::Impl::ElementType T;
    const T  One(1);
    copySum(x.left(), One, x.right(), y.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
// y = x1 - x2
//
template <typename VL, typename VR, typename VY>
typename RestrictTo<!ClosureType<OpSub, VL, VR>::isResidual,
         void>::Type
copy(const VectorClosure<OpSub, VL, VR> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(x, y);

    typedef typename VY::Impl::ElementType T;
    const T  MinusOne(-1);
    copySum(x.left(), MinusOne, x.right(), y.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  y = alpha*x
//
//  We evalute this with a scalSwitch
//  case 1: x is no closure
//  case 2: x is a closure

//
// case 1: x is no closure
//
template <typename ALPHA, typename VX, typename VY>
void
scalSwitch(const ALPHA &alpha, const Vector<VX> &x, Vector<VY> &y)
{
    using namespace DEBUGCLOSURE;

    if (identical(x.impl(), y.impl())) {
        scal(alpha, y.impl());
    } else {
//
//      Zero should have the same type as elements of y so that CBLAS
//      gets called if possible.
//
        typename Vector<VX>::Impl::ElementType  Zero(0);
        scal(Zero, y.impl());

        axpy(alpha, x.impl(), y.impl());
    }
}

//
// case 2: x is a closure
//
template <typename ALPHA, typename Op, typename L, typename R, typename VY>
void
scalSwitch(const ALPHA &alpha, const VectorClosure<Op, L, R> &x, Vector<VY> &y)
{
    copy(x, y.impl());
    scal(alpha, y.impl());
}

//
// entry point for switch
//
template <typename T, typename VX, typename VY>
void
copy(const VectorClosure<OpMult, ScalarValue<T>, VX> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(x, y);
    scalSwitch(x.left().value(), x.right(), y.impl());
    FLENS_BLASLOG_END;
}


//------------------------------------------------------------------------------
//
// y = x/alpha
//
// We evalute this with a rscalSwitch
// case 1: x is no closure
// case 2: x is a closure

//
// case 1: x is no closure
//
template <typename ALPHA, typename VX, typename VY>
void
rscalSwitch(const ALPHA &alpha, const Vector<VX> &x, Vector<VY> &y)
{
    using namespace DEBUGCLOSURE;

    if (identical(x.impl(), y.impl())) {
        rscal(alpha, y.impl());
    } else {
//
//      Zero should have the same type as elements of y so that CBLAS
//      gets called if possible.
//
        typename Vector<VX>::Impl::ElementType  Zero(0);
        scal(Zero, y.impl());

        raxpy(alpha, x.impl(), y.impl());
    }
}

//
// case 2: x is a closure
//
template <typename ALPHA, typename Op, typename L, typename R, typename VY>
void
rscalSwitch(const ALPHA &alpha, const VectorClosure<Op, L, R> &x, Vector<VY> &y)
{
    copy(x, y.impl());
    rscal(alpha, y.impl());
}

//
// entry point for switch
//
template <typename VX, typename T, typename VY>
void
copy(const VectorClosure<OpDiv, VX, ScalarValue<T> > &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(x, y);
    rscalSwitch(x.right().value(), x.left(), y.impl());
    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  y = conjugate(x)
//
template <typename VX, typename VY>
void
copy(const VectorClosureOpConj<VX> &, Vector<VY> &)
{
    ERROR_MSG("Lehn: Will be implemented on demand");
    ASSERT(0);
}

//------------------------------------------------------------------------------
//
//  y = beta*y + alpha*op(A)*x
//
template <typename VL, typename MVR, typename VY>
typename RestrictTo<ClosureType<OpAdd, VL, MVR>::isMatrixVectorProduct,
         void>::Type
copy(const VectorClosure<OpAdd, VL, MVR> &ypAx, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(ypAx, y);

    using namespace DEBUGCLOSURE;
//
//  check if y form rhs and lhs are identical
//
    typedef typename PruneScaling<VL>::Remainder    RVL;

    const RVL &_y = PruneScaling<VL>::getRemainder(ypAx.left());

    if (!identical(_y, y.impl())) {
        typedef typename VY::Impl::ElementType  TY;
        const TY  One(1);
        copySum(ypAx.left(), One, ypAx.right(), y.impl());
        FLENS_BLASLOG_END;
        return;
    }
//
//  get factor  beta
//
    typedef typename PruneScaling<VL>::ScalingType  SVL;
    const SVL &beta = PruneScaling<VL>::getFactor(ypAx.left());
//
//  Rest gets done by the mv switch
//
    typedef typename MVR::Left::ElementType TA;
    const auto &A = ypAx.right().left();
    const auto &x = ypAx.right().right();

    mvSwitch(NoTrans, TA(1), A, x, beta, y.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  y = A*x
//
template <typename ML, typename VR, typename VY>
typename RestrictTo<ClosureType<OpMult, ML, VR>::isMatrixVectorProduct,
         void>::Type
copy(const VectorClosure<OpMult, ML, VR> &Ax, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(Ax, y);

//
//  Call mv switch
//
    typedef typename ML::Impl::ElementType  TA;
    typedef typename VY::Impl::ElementType  TY;
    mvSwitch(NoTrans, TA(1), Ax.left(), Ax.right(), TY(0), y.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  y = x*A
//
template <typename VL, typename MR, typename VY>
typename RestrictTo<ClosureType<OpMult, VL, MR>::isVectorMatrixProduct,
         void>::Type
copy(const VectorClosure<OpMult, VL, MR> &xA, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(xA, y);

//
//  Call mv switch
//
    typedef typename MR::Impl::ElementType  TA;
    typedef typename VY::Impl::ElementType  TY;
    mvSwitch(Trans, TA(1), xA.right(), xA.left(), TY(0), y.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
// y = b - A*x
//
template <typename L, typename R, typename VY>
typename RestrictTo<ClosureType<OpSub, L, R>::isResidual,
         void>::Type
copy(const VectorClosure<OpSub, L, R> &bmAx, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(bmAx, y);

    const auto &b = bmAx.left();
    const auto &A = bmAx.right().left();
    const auto &x = bmAx.right().right();

    residual(b, A, x, y.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
// y = <Unknown Closure>
//
#ifndef FLENS_DEBUG_CLOSURES

template <typename Op, typename VL, typename VR, typename VY>
void
copy(const VectorClosure<Op, VL, VR> &, Vector<VY> &)
{
    ASSERT(0);
}

#else

template <typename Op, typename VL, typename VR, typename VY>
void
copy(const VectorClosure<Op, VL, VR> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_ERROR_COPY(x, y);
    ASSERT(0);
}

#endif

//== matrix closures ===========================================================
//
//  In the following comments op(X) denotes  X, X^T or X^H
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
    blas::copy(trans, A1, B);
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

//------------------------------------------------------------------------------
//
//  B = op(A1 + A2)
//
template <typename ML, typename MR, typename MB>
typename RestrictTo<!ClosureType<OpAdd, ML, MR>::isMatrixMatrixProduct,
         void>::Type
copy(Transpose trans, const MatrixClosure<OpAdd, ML, MR> &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_BEGIN_MCOPY(trans, A, B);

    typename MB::Impl::ElementType  One(1);
    copySum(trans, A.left(), One, A.right(), B.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  B = op(A1 - A2)
//
template <typename ML, typename MR, typename MB>
void
copy(Transpose trans, const MatrixClosure<OpSub, ML, MR> &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_BEGIN_MCOPY(trans, A, B);

    typename MB::Impl::ElementType  MinusOne(-1);
    copySum(trans, A.left(), MinusOne, A.right(), B.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  B = alpha*op(A)
//
//  We evalute this with a scalSwitch
//  case 1: A is no closure
//  case 2: A is a closure

//
// case 1: A is no closure
//
template <typename ALPHA, typename MA, typename MB>
void
scalSwitch(Transpose trans,
           const ALPHA &alpha, const Matrix<MA> &A, Matrix<MB> &B)
{
    using namespace DEBUGCLOSURE;

    if (identical(A.impl(), B.impl())) {
        scal(alpha, B.impl());
        if (trans!=NoTrans) {
            cotr(trans, B.impl());
        }
    } else {
//
//      Zero should have the same type as elements of y so that CBLAS
//      gets called if possible.
//
        typename Matrix<MA>::Impl::ElementType  Zero(0);
        scal(Zero, B.impl());

        axpy(trans, alpha, A.impl(), B.impl());
    }
}

//
// case 2: A is a closure
//
template <typename ALPHA, typename Op, typename L, typename R, typename MB>
void
scalSwitch(Transpose trans,
           const ALPHA &alpha, const MatrixClosure<Op, L, R> &A, Matrix<MB> &B)
{
    copy(trans, A, B.impl());
    scal(alpha, B.impl());
}

template <typename T, typename MA, typename MB>
void
copy(Transpose trans,
     const MatrixClosure<OpMult, ScalarValue<T>, MA> &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_BEGIN_COPY(A, B);
    scalSwitch(trans, A.left().value(), A.right(), B.impl());
    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  B = op(A)/alpha
//

//  We evalute this with a scalSwitch
//  case 1: A is no closure
//  case 2: A is a closure

//
// case 1: A is no closure
//
template <typename ALPHA, typename MA, typename MB>
void
rscalSwitch(Transpose trans,
            const ALPHA &alpha, const Matrix<MA> &A, Matrix<MB> &B)
{
    using namespace DEBUGCLOSURE;

    if (identical(A.impl(), B.impl())) {
        rscal(alpha, B.impl());
        if (trans!=NoTrans) {
            cotr(trans, B.impl());
        }
    } else {
//
//      Zero should have the same type as elements of y so that CBLAS
//      gets called if possible.
//
        typename Matrix<MA>::Impl::ElementType  Zero(0);
        scal(Zero, B.impl());

        raxpy(trans, alpha, A.impl(), B.impl());
    }
}

//
// case 2: A is a closure
//
template <typename ALPHA, typename Op, typename L, typename R, typename MB>
void
rscalSwitch(Transpose trans,
            const ALPHA &alpha, const MatrixClosure<Op, L, R> &A, Matrix<MB> &B)
{
    copy(trans, A, B.impl());
    rscal(alpha, B.impl());
}

template <typename MA, typename T, typename MB>
void
copy(Transpose trans,
     const MatrixClosure<OpDiv, MA, ScalarValue<T> > &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_BEGIN_COPY(A, B);
    rscalSwitch(trans, A.right().value(), A.left(), B.impl());
    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  B = op(A^T)
//
template <typename MA, typename MB>
void
copy(Transpose trans, const MatrixClosureOpTrans<MA> &A, Matrix<MB> &B)
{
    using namespace DEBUGCLOSURE;

    FLENS_BLASLOG_BEGIN_MCOPY(trans, A, B);

    Transpose _trans = Transpose(trans^Trans);

    copy(_trans, A.left(), B.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  C = beta*C + A*B
//
//  We evalute this with a mvSwitch
//  case 1: A is no closure
//  case 2: A is a scaling closure (i.e. scale*A)
//  case 3: A is some other closure
//  in all cases B is no closure.




//------------------------------------------------------------------------------
//
//  C = beta*C + A*B
//
template <typename ML, typename MR, typename MC>
typename RestrictTo<ClosureType<OpAdd, ML, MR>::isMatrixMatrixProduct,
         void>::Type
copy(Transpose trans, const MatrixClosure<OpAdd, ML, MR> &CPAB, Matrix<MC> &C)
{
    FLENS_BLASLOG_BEGIN_COPY(CPAB, C);

    using namespace DEBUGCLOSURE;
//
//  check if C form rhs and lhs are identical
//
    typedef typename PruneScaling<ML>::Remainder    RML;

    const RML &_C = PruneScaling<ML>::getRemainder(CPAB.left());

    if (trans!=NoTrans || !identical(_C, C.impl())) {
        typedef typename MC::Impl::ElementType  TC;
        const TC  One(1);
        copySum(trans, CPAB.left(), One, CPAB.right(), C.impl());
        FLENS_BLASLOG_END;
        return;
    }
//
//  get factor  beta
//
    typedef typename PruneScaling<ML>::ScalingType  SML;
    const SML &beta = PruneScaling<ML>::getFactor(CPAB.left());
//
//  Rest gets done by the mv switch
//
    typedef typename MR::Left::ElementType TA;
    const auto &A = CPAB.right().left();
    const auto &B = CPAB.right().right();

    if (trans==NoTrans || trans==Conj) {
        mmSwitch(trans, trans, TA(1), A, B, beta, C.impl());
    } else {
        mmSwitch(trans, trans, TA(1), B, A, beta, C.impl());
    }

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  C = A*B
//
template <typename MA, typename MB, typename MC>
typename RestrictTo<ClosureType<OpMult, MA, MB>::isMatrixMatrixProduct,
         void>::Type
copy(Transpose trans, const MatrixClosure<OpMult, MA, MB> &AB, Matrix<MC> &C)
{
    FLENS_BLASLOG_BEGIN_MCOPY(trans, AB, C);

//
//  Call mm switch
//
    typedef typename MA::Impl::ElementType  TA;
    typedef typename MC::Impl::ElementType  TC;

    const auto &A = AB.left();
    const auto &B = AB.right();

    if (trans==NoTrans || trans==Conj) {
        mmSwitch(trans, trans, TA(1), A, B, TC(0), C.impl());
    } else {
        mmSwitch(trans, trans, TA(1), B, A, TC(0), C.impl());
    }

    FLENS_BLASLOG_END;
}


//------------------------------------------------------------------------------
//
//  B = <Unknown Closure>
//
#ifndef FLENS_DEBUG_CLOSURES

template <typename Op, typename ML, typename MR, typename MB>
void
copy(Transpose, const MatrixClosure<Op, ML, MR> &, Matrix<MB> &)
{
    ERROR_MSG("B = <Unknown Closure>");
    ASSERT(0);
}

#else

template <typename Op, typename ML, typename MR, typename MB>
void
copy(Transpose trans, const MatrixClosure<Op, ML, MR> &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_ERROR_MCOPY(trans, A, B);
    ERROR_MSG("B = <Unknown Closure>");
    ASSERT(0);
}

#endif

//-- symmetric matrices --------------------------------------------------------
//
//  We just trans is NoTrans or Trans we simply ignore it
//
template <typename MA, typename MB>
void
copy(Transpose trans, const SymmetricMatrix<MA> &A, Matrix<MB> &B)
{
#   ifndef FLENS_DEBUG_CLOSURES
    ASSERT(trans==NoTrans || trans==Trans);
#   else
    if (trans!=NoTrans && trans!=Trans) {
        typedef typename MA::ElementType TA;

        GeMatrix<FullStorage<TA> >  _A = A;
        FLENS_BLASLOG_TMP_ADD(_A);

        copy(trans, _A, B.impl());

        FLENS_BLASLOG_TMP_REMOVE(_A, A);
        return;
    }
#   endif
    copy(A.impl(), B.impl());
}

template <typename MA, typename MB>
typename RestrictTo<!IsSymmetricMatrix<MA>::value,
         void>::Type
copy(Transpose trans, const Matrix<MA> &A, SymmetricMatrix<MB> &B)
{
    copy(trans, A.impl(), B.impl());
}

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_COPY_TCC
