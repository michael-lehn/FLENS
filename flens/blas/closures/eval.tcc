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

#ifndef FLENS_BLAS_CLOSURES_EVAL_TCC
#define FLENS_BLAS_CLOSURES_EVAL_TCC 1

#include <flens/blas/closures/debugclosure.h>
#include <flens/blas/closures/prunematrixclosure.h>
#include <flens/blas/closures/prunevectorclosure.h>
#include <flens/blas/debugmacro.h>
#include <flens/blas/level1/level1.h>
#include <flens/blas/level2/level2.h>

//-- blas entry points----------------------------------------------------------

namespace flens {

//-- vector closures
template <typename VX, typename VY>
void
copy(const Vector<VX> &x, Vector<VY> &y)
{
    FLENS_CLOSURELOG_ADD_ENTRY_ASSIGNMENT(x, y);

    blas::copy(x.impl(), y.impl());

    FLENS_CLOSURELOG_END_ENTRY;
}

template <typename ALPHA, typename VX, typename VY>
void
axpy(const ALPHA &alpha, const Vector<VX> &x, Vector<VY> &y)
{
    FLENS_CLOSURELOG_ADD_ENTRY_PLUSASSIGNMENT(alpha, x, y);

    blas::axpy(alpha, x.impl(), y.impl());

    FLENS_CLOSURELOG_END_ENTRY;
}

//-- matrix closures
template <typename MA, typename MB>
void
copy(cxxblas::Transpose trans,
     const Matrix<MA> &A, Matrix<MB> &B)
{
    // TODO: take trans into account
    FLENS_CLOSURELOG_ADD_ENTRY_ASSIGNMENT(A, B);

    blas::copy(trans, A.impl(), B.impl());

    FLENS_CLOSURELOG_END_ENTRY;
}

template <typename ALPHA, typename MA, typename MB>
void
axpy(cxxblas::Transpose trans,
     const ALPHA &alpha, const Matrix<MA> &A, Matrix<MB> &B)
{
    // TODO: take trans into account
    FLENS_CLOSURELOG_ADD_ENTRY_PLUSASSIGNMENT(alpha, A, B);

    blas::axpy(trans, alpha, A.impl(), B.impl());

    FLENS_CLOSURELOG_END_ENTRY;
}

} // namespace flens

namespace flens { namespace blas {

//-- vector closures -----------------------------------------------------------
//-- copy
// y = alpha*x
// (eval-switch case 1: x is a closure)
template <typename ALPHA, typename Op, typename L, typename R, typename VY>
void
copyScal(const ALPHA &alpha, const VectorClosure<Op, L, R> &x, Vector<VY> &y)
{
    copy(x, y.impl());
    scal(alpha, y.impl());
}

// y = alpha*x
// (eval-switch case 2: x is no closure)
template <typename ALPHA, typename VX, typename VY>
void
copyScal(const ALPHA &alpha, const Vector<VX> &x, Vector<VY> &y)
{
    if (ADDRESS(x)==ADDRESS(y)) {
        scal(alpha, y.impl());
    } else {
        scal(ALPHA(0), y.impl());
        axpy(alpha, x.impl(), y.impl());
    }
}

// y = alpha*x
// (entry point for eval-switch)
template <typename T, typename VX, typename VY>
void
copy(const VectorClosure<OpMult, ScalarValue<T>, VX> &ax, Vector<VY> &y)
{
    /*
    if (ADDRESS(ax.right())==ADDRESS(y)) {
        scal(ax.left().value(), y.impl());
    } else {
        scal(T(0), y.impl());
        axpy(ax.left().value(), ax.right(), y.impl());
    }
    */
    copyScal(ax.left().value(), ax.right(), y.impl());
}

// y = x1 + x2
template <typename VL, typename VR, typename VY>
void
copy(const VectorClosure<OpAdd, VL, VR> &x, Vector<VY> &y)
{
    FLENS_CLOSURELOG_ADD_ENTRY_COPY(x, y);

    ASSERT(!DebugClosure::search(x.right(), ADDRESS(y)));
    copy(x.left(), y.impl());

    typedef typename VY::Impl::ElementType TY;
    axpy(TY(1), x.right(), y.impl());

    FLENS_CLOSURELOG_END_ENTRY;
}

// y = x1 - x2
template <typename VL, typename VR, typename VY>
void
copy(const VectorClosure<OpSub, VL, VR> &x, Vector<VY> &y)
{
    ASSERT(!DebugClosure::search(x.right(), ADDRESS(y)));

    copy(x.left(), y.impl());
    typedef typename VY::Impl::ElementType TY;
    axpy(TY(-1), x.right(), y.impl());
}

// y = A*x
template <typename MA, typename VX, typename VY>
void
copy(const VectorClosure<OpMult, MA, VX> &Ax, Vector<VY> &y)
{
    typedef typename VY::Impl::ElementType TY;
    mv(cxxblas::NoTrans, TY(1), Ax.left(), Ax.right(), TY(0), y.impl());
}

//-- axpy
// y += alpha*x
template <typename ALPHA, typename T, typename VX, typename VY>
void
axpy(const ALPHA &alpha,
     const VectorClosure<OpMult, ScalarValue<T>, VX> &ax, Vector<VY> &y)
{
    axpy(alpha*ax.left().value(), ax.right(), y.impl());
}

// y += x1 + x2
template <typename ALPHA, typename VL, typename VR, typename VY>
void
axpy(const ALPHA &alpha,
     const VectorClosure<OpAdd, VL, VR> &x, Vector<VY> &y)
{
    ASSERT(!DebugClosure::search(x.right(), ADDRESS(y)));
    FLENS_CLOSURELOG_ADD_ENTRY_AXPY(alpha, x, y);

    axpy(alpha, x.left(), y.impl());
    axpy(alpha, x.right(), y.impl());

    FLENS_CLOSURELOG_END_ENTRY;
}

// y += x1 - x2
template <typename ALPHA, typename VL, typename VR, typename VY>
void
axpy(const ALPHA &alpha,
     const VectorClosure<OpSub, VL, VR> &x, Vector<VY> &y)
{
    ASSERT(!DebugClosure::search(x.right(), ADDRESS(y)));

    axpy(alpha, x.left(), y.impl());
    axpy(-alpha, x.right(), y.impl());
}

// y += A*x
template <typename ALPHA, typename MA, typename VX, typename VY>
void
axpy(const ALPHA &alpha,
     const VectorClosure<OpMult, MA, VX> &Ax, Vector<VY> &y)
{
    typedef typename VY::Impl::ElementType TY;
    mv(cxxblas::NoTrans, alpha, Ax.left(), Ax.right(), TY(1), y.impl());
}

//-- mv
// y = x*A  ->  y = A'*x
template <typename ALPHA, typename VX, typename MA, typename BETA, typename VY>
void
mv(cxxblas::Transpose trans,
   const ALPHA &alpha, const Vector<VX> &x, const Matrix<MA> &A,
   const BETA &beta, Vector<VY> &y)
{
    trans = cxxblas::Transpose(cxxblas::Trans^trans);
    mv(trans, alpha, A.impl(), x.impl(), beta, y.impl());
}

template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
void
mv(cxxblas::Transpose trans,
   const ALPHA &alpha, const Matrix<MA> &A, const Vector<VX> &x,
   const BETA &beta, Vector<VY> &y)
{
    typedef PruneMatrixClosure<typename MA::Impl> PMC;
    typedef PruneVectorClosure<typename VX::Impl> PVC;
    trans = PVC::updateTranspose(trans);
    trans = PMC::updateTranspose(trans);
    mv(trans, 
       PMC::updateScalingFactor(PVC::updateScalingFactor(alpha, x.impl()),
                                A.impl()),
       PMC::remainder(A.impl()),
       PVC::remainder(x.impl()),
       beta,
       y.impl());
}

template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
void
mv(cxxblas::Transpose trans,
   const ALPHA &alpha, const HermitianMatrix<MA> &A, const Vector<VX> &x,
   const BETA &beta, Vector<VY> &y)
{
    ASSERT(trans==cxxblas::NoTrans);
    mv(alpha, A.impl(), x.impl(), beta, y.impl());
}

template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
void
mv(cxxblas::Transpose trans,
   const ALPHA &alpha, const SymmetricMatrix<MA> &A, const Vector<VX> &x,
   const BETA &beta, Vector<VY> &y)
{
    ASSERT(trans==cxxblas::NoTrans);
    mv(alpha, A.impl(), x.impl(), beta, y.impl());
}

template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
void
mv(cxxblas::Transpose trans,
   const ALPHA &alpha, const TriangularMatrix<MA> &A, const Vector<VX> &x,
   const BETA &beta, Vector<VY> &y)
{
    ASSERT(alpha==ALPHA(1));
    ASSERT(beta==BETA(0));
    ASSERT(ADDRESS(x)==ADDRESS(y));
    mv(trans, A.impl(), y.impl());
}

//-- matrix closures -----------------------------------------------------------
//-- copy
// B = op(A)   (error-handle for unkown closure)
template <typename Op, typename MA1, typename MA2, typename MB>
void
copy(cxxblas::Transpose trans,
     const MatrixClosure<Op, MA1, MA2> &A, Matrix<MB> &B)
{
    std::cerr << "ERROR: unkown Op-type" << std::endl;
    ASSERT(0);
}

// B = A'
template <typename MA, typename MB>
void
copy(cxxblas::Transpose trans,
     const MatrixClosure<OpTrans, MA, MA> &At, Matrix<MB> &B)
{
    trans = cxxblas::Transpose(trans^cxxblas::Trans);
    if (ADDRESS(At.left())==ADDRESS(B)) {
        ASSERT(trans==cxxblas::NoTrans);
    } else {
        copy(trans, At.left(), B.impl());
    }
}

// B = alpha*A
template <typename T, typename MA, typename MB>
void
copy(cxxblas::Transpose trans,
     const MatrixClosure<OpMult, ScalarValue<T>, MA> &aA, Matrix<MB> &B)
{
    if (ADDRESS(aA.right())==ADDRESS(B)) {
        ASSERT(trans==cxxblas::NoTrans);
        scal(aA.left().value(), B.impl());
    } else {
        scal(T(0), B.impl());
        axpy(trans, aA.left().value(), aA.right(), B.impl());
    }
}

// B = A1 + A2
template <typename ML, typename MR, typename MB>
void
copy(cxxblas::Transpose trans,
     const MatrixClosure<OpAdd, ML, MR> &A, Matrix<MB> &B)
{
    FLENS_CLOSURELOG_ADD_ENTRY_COPY(A, B);

    ASSERT(!DebugClosure::search(A.right(), ADDRESS(B)));
    copy(trans, A.left(), B.impl());

    typedef typename MB::Impl::ElementType TB;
    axpy(trans, TB(1), A.right(), B.impl());

    FLENS_CLOSURELOG_END_ENTRY;
}

//-- axpy
// B += op(A)   (error-handle for unkown closure)
template <typename ALPHA, typename Op, typename MA1, typename MA2, typename MB>
void
axpy(cxxblas::Transpose trans, const ALPHA &alpha,
     const MatrixClosure<Op, MA1, MA2> &A, Matrix<MB> &B)
{
    std::cerr << "ERROR: unkown Op-type" << std::endl;
    ASSERT(0);
}

// B += A'
template <typename ALPHA, typename MA, typename MB>
void
axpy(cxxblas::Transpose trans, const ALPHA &alpha,
     const MatrixClosure<OpTrans, MA, MA> &At, Matrix<MB> &B)
{
    trans = cxxblas::Transpose(trans^cxxblas::Trans);
    axpy(trans, alpha, At.left(), B.impl());
}

// B += alpha*A
template <typename ALPHA, typename T, typename MA, typename MB>
void
axpy(cxxblas::Transpose trans, const ALPHA &alpha,
     const MatrixClosure<OpMult, ScalarValue<T>, MA> &aA, Matrix<MB> &B)
{
    axpy(trans, alpha*aA.left().value(), aA.right(), B.impl());
}

// B += A1 + A2
template <typename ALPHA, typename ML, typename MR, typename MB>
void
axpy(cxxblas::Transpose trans, const ALPHA &alpha,
     const MatrixClosure<OpAdd, ML, MR> &A, Matrix<MB> &B)
{
    ASSERT(!DebugClosure::search(A.right(), ADDRESS(B)));
    FLENS_CLOSURELOG_ADD_ENTRY_AXPY(alpha, A, B);

    axpy(trans, alpha, A.left(), B.impl());
    axpy(trans, alpha, A.right(), B.impl());

    FLENS_CLOSURELOG_END_ENTRY;
}

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_EVAL_TCC
