/*
 *   Copyright (c) 2012, Michael Lehn, Klaus Pototzky
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

#ifndef FLENS_BLAS_LEVEL1_AXPY_TCC
#define FLENS_BLAS_LEVEL1_AXPY_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/closures/closures.h>
#include <flens/blas/level1/level1.h>
#include <flens/typedefs.h>
#include <ulmblas/cxxblas.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//-- axpy
template <typename ALPHA, typename VX, typename VY>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
axpy(const ALPHA &alpha, const VX &x, VY &&y)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_AXPY(alpha, x, y);

    if (y.length()==0) {
//
//      We allow  y += alpha*x  for an empty vector y
//
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;
        const T  Zero(0);

        y.resize(x, Zero);
    }
    ASSERT(y.length()==x.length());

    cxxblas::axpy(x.length(), alpha,
                  x.data(), x.stride(),
                  y.data(), y.stride());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- acxpy
template <typename ALPHA, typename VX, typename VY>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
acxpy(const ALPHA &alpha, const VX &x, VY &&y)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_AXPY(alpha, x, y);

    if (y.length()==0) {
//
//      We allow  y += alpha*conjugate(x)  for an empty vector y
//
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;
        const T  Zero(0);

        y.resize(x, Zero);
    }
    ASSERT(y.length()==x.length());

    cxxblas::acxpy(x.length(), alpha,
                   x.data(), x.stride(),
                   y.data(), y.stride());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- geaxpy
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
axpy(Transpose transposeA, const ALPHA &alpha, const MA &A, MB &&B)
{
    typedef typename RemoveRef<MB>::Type   MatrixB;

    const bool transA = (transposeA==Trans || transposeA==ConjTrans);
    const bool conjA  = (transposeA==Conj || transposeA==ConjTrans);

    if (B.numRows()==0 || B.numCols()==0) {
//
//      So we allow  B += alpha*A  for an empty matrix B
//
        typedef typename MatrixB::ElementType  T;
        const T  Zero(0);

        if (!transA) {
            B.resize(A.numRows(), A.numCols(),
                     A.firstRow(), A.firstCol(),
                     Zero);
        } else {
            B.resize(A.numCols(), A.numRows(),
                     A.firstCol(), A.firstRow(),
                     Zero);
        }
    }

#   ifndef NDEBUG
    if (!transA) {
        ASSERT((A.numRows()==B.numRows()) && (A.numCols()==B.numCols()));
    } else {
        ASSERT((A.numRows()==B.numCols()) && (A.numCols()==B.numRows()));
    }
#   endif

#   ifndef FLENS_DEBUG_CLOSURES
#   ifndef NDEBUG
    if (transA) {
        ASSERT(!DEBUGCLOSURE::identical(A, B));
    }
#   endif
#   else
//
//  If A and B are identical a temporary is needed if we want to use axpy
//  for B += alpha*A^T or B+= alpha*A^H
//
    if (transA && DEBUGCLOSURE::identical(A, B)) {
        typename Result<MA>::Type A_ = A;
        FLENS_BLASLOG_TMP_ADD(A_);

        copy(transposeA, A, A_);
        axpy(NoTrans, alpha, A, B);

        FLENS_BLASLOG_TMP_REMOVE(A_, A);
        return;
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MAXPY(transposeA, alpha, A, B);

    cxxblas::geaxpy(B.numRows(), B.numCols(),
                    alpha,
                    transA, conjA, A.data(), A.strideRow(), A.strideCol(),
                    B.data(), B.strideRow(), B.strideCol());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_AXPY_TCC
