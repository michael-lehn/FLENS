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

#ifndef FLENS_BLAS_LEVEL1_RAXPY_TCC
#define FLENS_BLAS_LEVEL1_RAXPY_TCC 1

#include <cxxblas/cxxblas.h>
#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/closures/closures.h>
#include <flens/blas/level1/level1.h>
#include <flens/typedefs.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//-- raxpy
template <typename ALPHA, typename VX, typename VY>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
raxpy(const ALPHA &alpha, const VX &x, VY &&y)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_RAXPY(alpha, x, y);

    if (y.length()==0) {
//
//      So we allow  y += 1/alpha*x  for an empty vector y
//
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;
        const T  Zero(0);

        y.resize(x, Zero);
    }
    ASSERT(y.length()==x.length());

#   ifdef HAVE_CXXBLAS_RAXPY
    cxxblas::raxpy(x.length(), alpha,
                   x.data(), x.stride(),
                   y.data(), y.stride());
#   else
    ASSERT(0);
#   endif
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- tinyvector raxpy
//
template <typename ALPHA, typename VX, typename VY>
typename RestrictTo<IsTinyVector<VX>::value
                 && IsTinyVector<VY>::value,
         void>::Type
raxpy(const ALPHA &alpha, const VX &x, VY &&y)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_RAXPY(alpha, x, y);

    typedef typename VX::ElementType      TX;
    typedef typename RemoveRef<VY>::Type  VectorY;
    typedef typename VectorY::ElementType TY;

    const int lengthX = VX::Engine::length;
    const int strideX = VX::Engine::stride;
    const int lengthY = VectorY::Engine::length;
    const int strideY = VectorY::Engine::stride;

    ASSERT(lengthY == lengthX && lengthY != 0);

#   ifdef HAVE_CXXBLAS_RAXPY
    cxxblas::raxpy(lengthX, alpha, x.data(), strideX, y.data(), strideY);
#   else
    ASSERT(0);
#   endif
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- geraxpy
//
//  B += A/alpha
//
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
raxpy(Transpose trans, const ALPHA &alpha, const MA &A, MB &&B)
{
    typedef typename RemoveRef<MB>::Type   MatrixB;

    if (B.numRows()==0 || B.numCols()==0) {
//
//      So we allow  B += 1/alpha*A  for an empty matrix B
//
        typedef typename MatrixB::ElementType  T;
        const T  Zero(0);

        if ((trans==NoTrans) || (trans==Conj)) {
            B.resize(A.numRows(), A.numCols(), Zero);
        } else {
            B.resize(A.numCols(), A.numRows(), Zero);
        }
    }

#   ifndef NDEBUG
    if ((trans==NoTrans) || (trans==Conj)) {
        ASSERT((A.numRows()==B.numRows()) && (A.numCols()==B.numCols()));
    } else {
        ASSERT((A.numRows()==B.numCols()) && (A.numCols()==B.numRows()));
    }
#   endif


    trans = (A.order()==B.order())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);


#   ifndef FLENS_DEBUG_CLOSURES
#   ifndef NDEBUG
    if (trans==Trans || trans==ConjTrans) {
        ASSERT(!DEBUGCLOSURE::identical(A, B));
    }
#   endif
#   else
    typedef typename RemoveRef<MA>::Type  MatrixA;
//
//  If A and B are identical a temporary is needed if we want to use axpy
//  for B += A^T/alpha or B+= A^H/alpha
//
    if ((trans==Trans || trans==ConjTrans) && DEBUGCLOSURE::identical(A, B)) {

        typename Result<MatrixA>::Type A_ = A;
        FLENS_BLASLOG_TMP_ADD(A_);

        axpy(trans, alpha, A_, B);

        FLENS_BLASLOG_TMP_REMOVE(A_, A);
        return;
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MRAXPY(trans, alpha, A, B);

#   ifdef HAVE_CXXBLAS_GERAXPY
    geraxpy(B.order(), trans, B.numRows(), B.numCols(), alpha,
            A.data(), A.leadingDimension(),
            B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}


} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_RAXPY_TCC
