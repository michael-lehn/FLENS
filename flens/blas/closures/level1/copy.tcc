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

#ifndef FLENS_BLAS_CLOSURES_LEVEL1_COPY_TCC
#define FLENS_BLAS_CLOSURES_LEVEL1_COPY_TCC 1

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

//------------------------------------------------------------------------------
//
// y = x1 + x2
//
template <typename VL, typename VR, typename VY>
typename RestrictTo<VCDefaultEval<OpAdd, VL, VR>::value
                 && IsVector<VL>::value
                 && IsVector<VR>::value,
         void>::Type
copy(const VectorClosure<OpAdd, VL, VR> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(false, x, y);

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
typename RestrictTo<VCDefaultEval<OpSub, VL, VR>::value
                 && IsVector<VL>::value
                 && IsVector<VR>::value,
         void>::Type
copy(const VectorClosure<OpSub, VL, VR> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(false, x, y);

    typedef typename VY::Impl::ElementType T;
    const T  MinusOne(-1);
    copySum(x.left(), MinusOne, x.right(), y.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  y = alpha*x
//

template <typename SV, typename VX, typename VY>
typename RestrictTo<VCDefaultEval<OpMult, SV, VX>::value
                 && IsScalarValue<SV>::value
                 && IsVector<VX>::value,
         void>::Type
copy(const VectorClosure<OpMult, SV, VX> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(false, x, y);
    copyScal(x.left().value(), x.right(), y.impl());
    FLENS_BLASLOG_END;
}


//------------------------------------------------------------------------------
//
// y = x/alpha
//

//
// entry point for switch
//
template <typename VX, typename SV, typename VY>
typename RestrictTo<VCDefaultEval<OpDiv, VX, SV>::value
                 && IsVector<VX>::value
                 && IsScalarValue<SV>::value,
         void>::Type
copy(const VectorClosure<OpDiv, VX, SV> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(false, x, y);
    copyRScal(x.right().value(), x.left(), y.impl());
    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  y = conjugate(x)
//
template <typename VX, typename VY>
typename RestrictTo<DefaultEval<VectorClosureOpConj<VX> >::value
                 && IsVector<VX>::value,
         void>::Type
copy(const VectorClosureOpConj<VX> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(false, x, y);
    copyConj(x.left(), y.impl());
    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  y = A*x
//
template <typename ML, typename VR, typename VY>
typename RestrictTo<VCDefaultEval<OpMult, ML, VR>::value
                 && IsMatrix<ML>::value
                 && IsVector<VR>::value,
         void>::Type
copy(const VectorClosure<OpMult, ML, VR> &Ax, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(false, Ax, y);

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
typename RestrictTo<VCDefaultEval<OpMult, VL, MR>::value
                 && IsVector<VL>::value
                 && IsMatrix<MR>::value,
         void>::Type
copy(const VectorClosure<OpMult, VL, MR> &xA, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(false, xA, y);

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
// y = <Unknown Closure>
//

template <typename Op, typename VL, typename VR, typename VY>
void
copy(const VectorClosure<Op, VL, VR> &FLENS_BLASLOG_VARDECL(x),
     Vector<VY>                      &FLENS_BLASLOG_VARDECL(y))
{
    FLENS_BLASLOG_ERROR_COPY(x, y);

    ASSERT(0);
}

//== matrix closures ===========================================================
//
//  In the following comments op(X) denotes  X, X^T or X^H
//

//------------------------------------------------------------------------------
//
//  B = op(A1 + A2)
//
template <typename ML, typename MR, typename MB>
typename RestrictTo<MCDefaultEval<OpAdd, ML, MR>::value
                 && IsMatrix<ML>::value
                 && IsMatrix<MR>::value,
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
typename RestrictTo<MCDefaultEval<OpSub, ML, MR>::value
                 && IsMatrix<ML>::value
                 && IsMatrix<MR>::value,
         void>::Type
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

template <typename SV, typename MA, typename MB>
typename RestrictTo<MCDefaultEval<OpMult, SV, MA>::value
                 && IsScalarValue<SV>::value
                 && IsMatrix<MA>::value,
         void>::Type
copy(Transpose trans, const MatrixClosure<OpMult, SV, MA> &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_BEGIN_COPY(false, A, B);
    copyScal(trans, A.left().value(), A.right(), B.impl());
    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  B = op(A)/alpha
//

template <typename MA, typename SV, typename MB>
typename RestrictTo<MCDefaultEval<OpDiv, MA, SV>::value
                 && IsMatrix<MA>::value
                 && IsScalarValue<SV>::value,
         void>::Type
copy(Transpose trans, const MatrixClosure<OpDiv, MA, SV> &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_BEGIN_COPY(false, A, B);
    copyRScal(trans, A.right().value(), A.left(), B.impl());
    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  B = op(conjugate(A))
//
template <typename MA, typename MB>
typename RestrictTo<DefaultEval<MatrixClosureOpConj<MA> >::value
                 && IsMatrix<MA>::value,
         void>::Type
copy(Transpose trans, const MatrixClosureOpConj<MA> &A, Matrix<MB> &B)
{
    using namespace DEBUGCLOSURE;

    FLENS_BLASLOG_BEGIN_MCOPY(trans, A, B);

    Transpose trans_ = Transpose(trans^Conj);

    copy(trans_, A.left(), B.impl());

    FLENS_BLASLOG_END;
}

//------------------------------------------------------------------------------
//
//  B = op(A^T)
//
template <typename MA, typename MB>
typename RestrictTo<DefaultEval<MatrixClosureOpTrans<MA> >::value
                 && IsMatrix<MA>::value,
         void>::Type
copy(Transpose trans, const MatrixClosureOpTrans<MA> &A, Matrix<MB> &B)
{
    using namespace DEBUGCLOSURE;

    FLENS_BLASLOG_BEGIN_MCOPY(trans, A, B);

    Transpose trans_ = Transpose(trans^Trans);

    copy(trans_, A.left(), B.impl());

    FLENS_BLASLOG_END;
}


//------------------------------------------------------------------------------
//
//  C = A*B
//
template <typename MA, typename MB, typename MC>
typename RestrictTo<MCDefaultEval<OpMult, MA, MB>::value
                 && IsMatrix<MA>::value
                 && IsMatrix<MB>::value,
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
#ifdef FLENS_DEBUG_CLOSURES

template <typename Op, typename ML, typename MR, typename MB>
void
copy(Transpose DEBUG_VAR(trans), const MatrixClosure<Op, ML, MR> &DEBUG_VAR(A),
     Matrix<MB> &DEBUG_VAR(B) )
{
    FLENS_BLASLOG_ERROR_MCOPY(trans, A, B);
    ERROR_MSG("B = <Unknown Closure>");

    ASSERT(0);
}

#endif

//-- hermitian matrices --------------------------------------------------------

template <typename MA, typename MB>
typename RestrictTo<IsHermitianMatrix<MA>::value,
         void>::Type
copy(Transpose DEBUG_VAR(trans), const MA &A, Matrix<MB> &B)
{
#   ifndef FLENS_DEBUG_CLOSURES
    ASSERT(trans==NoTrans || trans==ConjTrans);
#   else
    if (trans!=NoTrans && trans!=Trans) {
        typedef typename MA::ElementType TA;

        GeMatrix<FullStorage<TA> >  A_ = A;
        FLENS_BLASLOG_TMP_ADD(A_);

        copy(trans, A_, B.impl());

        FLENS_BLASLOG_TMP_REMOVE(A_, A);
        return;
    }
#   endif
    copy(A.impl(), B.impl());
}

//-- symmetric matrices --------------------------------------------------------

template <typename MA, typename MB>
typename RestrictTo<IsSymmetricMatrix<MA>::value,
         void>::Type
copy(Transpose DEBUG_VAR(trans), const MA &A, Matrix<MB> &B)
{
#   ifndef FLENS_DEBUG_CLOSURES
    ASSERT(trans==NoTrans || trans==Trans);
#   else
    if (trans!=NoTrans && trans!=Trans) {
        typedef typename MA::ElementType TA;

        GeMatrix<FullStorage<TA> >  A_ = A;
        FLENS_BLASLOG_TMP_ADD(A_);

        copy(trans, A_, B.impl());

        FLENS_BLASLOG_TMP_REMOVE(A_, A);
        return;
    }
#   endif
    copy(A.impl(), B.impl());
}

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_LEVEL1_COPY_TCC
