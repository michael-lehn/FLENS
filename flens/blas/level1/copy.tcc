/*
 *   Copyright (c) 2009, Michael Lehn
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

#ifndef FLENS_BLAS_LEVEL1_COPY_TCC
#define FLENS_BLAS_LEVEL1_COPY_TCC 1

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

//-- copy
template <typename VX, typename VY>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
copy(const VX &x, VY &&y)
{
    copy(false, x, y);
}

template <typename VX, typename VY>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
copy(bool conjX, const VX &x, VY &&y)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_COPY(conjX, x, y);

//
//  Resize left hand size if needed.  This is *usually* only allowed
//  when the left hand side is an empty vector (such that it is no actual
//  resizing but rather an initialization).
//
    if (y.length()!=x.length()) {
#       ifndef FLENS_DEBUG_CLOSURES
        ASSERT(y.length()==0);
#       else
        if (y.length()!=0) {
            FLENS_BLASLOG_RESIZE_VECTOR(y, x.length());
        }
#       endif
        y.resize(x);
    }

    cxxblas::copy(x.length(),
                  conjX,  x.data(), x.stride(),
                  y.data(), y.stride());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- gecopy
template <typename MA, typename MB>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
copy(Transpose transposeA, const MA &A, MB &&B)
{
    const bool transA = (transposeA==Trans || transposeA==ConjTrans);
    const bool conjA  = (transposeA==Conj || transposeA==ConjTrans);

//
//  check if this is an inplace transpose of A
//
    if (transA) {
        if (DEBUGCLOSURE::identical(A, B)) {
#           ifndef FLENS_DEBUG_CLOSURES
//
//          temporaries are not allowed
//
            ASSERT(B.numRows()==B.numCols());
            cxxblas::gecotr(B.numRows(), transA, conjA,
                            B.data(), B.strideRow(), B.strideCol());
            return;
#           else
//
//          temporaries are allowed: check if this requires a temporary
//
            if (A.numRows()!=A.numCols()) {
                typedef typename RemoveRef<MA>::Type   MatrixA;

                typename Result<MatrixA>::Type A_ = A;
                FLENS_BLASLOG_TMP_ADD(A_);

                copy(transposeA, A_, B);

                FLENS_BLASLOG_TMP_REMOVE(A_, A);
                return;
            } else {
//
//              otherwise perform inplace transpose
//
                cxxblas::gecotr(B.numRows(), transA, conjA,
                                B.data(), B.strideRow(), B.strideCol());
                return;
            }
#           endif
        }
    }

//
//  Resize left hand size if needed.  This is *usually* only allowed
//  when the left hand side is an empty matrix (such that it is no actual
//  resizing but rather an initialization).
//
    if (!transA) {
        if ((A.numRows()!=B.numRows()) || (A.numCols()!=B.numCols())) {
#           ifndef FLENS_DEBUG_CLOSURES
            ASSERT(B.numRows()==0 || B.numCols()==0);
#           else
            if (B.numRows()!=0 || B.numCols()!=0) {
                FLENS_BLASLOG_RESIZE_MATRIX(B, A.numRows(), A.numCols());
            }
#           endif
            B.resize(A);
        }
    } else {
        if ((A.numRows()!=B.numCols())  || (A.numCols()!=B.numRows())) {
#           ifndef FLENS_DEBUG_CLOSURES
            ASSERT(B.numRows()==0 || B.numCols()==0);
#           else
            if (B.numRows()!=0 || B.numCols()!=0) {
                FLENS_BLASLOG_RESIZE_MATRIX(B, A.numCols(), A.numRows());
            }
#           endif
            B.resize(A.numCols(), A.numRows(),
                     A.firstCol(), A.firstRow());
        }
    }

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MCOPY(transposeA, A, B);

    cxxblas::gecopy(B.numRows(), B.numCols(),
                    transA, conjA, A.data(), A.strideRow(), A.strideCol(),
                    B.data(), B.strideRow(), B.strideCol());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- trcopy
template <typename MA, typename MB>
typename RestrictTo<IsTrMatrix<MA>::value
                 && IsTrMatrix<MB>::value,
         void>::Type
copy(Transpose transposeA, const MA &A, MB &&B)
{
    const bool transA    = (transposeA==Trans || transposeA==ConjTrans);
    const bool conjA     = (transposeA==Conj || transposeA==ConjTrans);
    const bool lowerA    = (A.upLo()==Lower);
    const bool unitDiagA = (A.diag()==Unit);

    ASSERT(A.diag()==B.diag());
    ASSERT((!transA && A.upLo()==B.upLo()) || (transA && A.upLo()!=B.upLo()));

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MCOPY(transposeA, A, B);

    cxxblas::trcopy(B.numRows(), B.numCols(),
                    lowerA, transA, conjA, unitDiagA,
                    A.data(), A.strideRow(), A.strideCol(),
                    B.data(), B.strideRow(), B.strideCol());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_COPY_TCC
