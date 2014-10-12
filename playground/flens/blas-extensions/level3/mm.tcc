/*
 *   Copyright (c) 2013, Michael Lehn, Klaus Pototzky
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

#ifndef PLAYGROUND_FLENS_BLASEXTENSIONS_LEVEL3_MM_TCC
#define PLAYGROUND_FLENS_BLASEXTENSIONS_LEVEL3_MM_TCC

#include <flens/blas/closures/closures.h>
#include <flens/blas/level3/level3.h>
#include <flens/typedefs.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas { namespace extensions {

//== GeneralMatrix - GeneralMatrix products ====================================

//-- gemm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsGeMatrix<MB>::value
                 && IsGeMatrix<MC>::value,
         void>::Type
mm_strassen(Transpose        transposeA,
            Transpose        transposeB,
            const ALPHA      &alpha,
            const MA         &A,
            const MB         &B,
            const BETA       &beta,
            MC               &&C)
{
    typedef typename RemoveRef<MC>::Type MatrixC;
    typedef typename MatrixC::IndexType  IndexType;

    const bool noTransA = (transposeA==NoTrans || transposeA==Conj);
    const bool noTransB = (transposeB==NoTrans || transposeB==Conj);

#   ifndef NDEBUG
    int kA = (noTransA) ? A.numCols() : A.numRows();
    int kB = (noTransB) ? B.numRows() : B.numCols();
    ASSERT(kA==kB);
#   endif

    IndexType m = (noTransA) ? A.numRows() : A.numCols();
    IndexType n = (noTransB) ? B.numCols() : B.numRows();
    IndexType k = (noTransA) ? A.numCols() : A.numRows();

    if (C.order()!=A.order()) {
        transposeA = Transpose(transposeA ^ Trans);
    }
    if (C.order()!=B.order()) {
        transposeB = Transpose(transposeB ^ Trans);
    }

#   ifndef NDEBUG
    if (beta!=BETA(0)) {
        if (C.numRows()!=0 && C.numCols()!=0) {
            ASSERT(C.numRows()==m && C.numCols()==n);
        }
    }
#   endif

    if ((C.numRows()!=m) || (C.numCols()!=n)) {
        ASSERT(C.numRows()==0 && C.numCols()==0);
        C.resize(m, n);
    }

#   ifndef FLENS_DEBUG_CLOSURES
    ASSERT(!DEBUGCLOSURE::identical(A, C));
    ASSERT(!DEBUGCLOSURE::identical(B, C));
#   else
//
//  If A or B is identical with C we copy C into a temporary first.  Then
//  we compute the matrix-matrix product and afterwards copy the result into C.
//
    if (DEBUGCLOSURE::identical(A, C) || DEBUGCLOSURE::identical(B, C)) {
        typename MatrixC::NoView C_;
        FLENS_BLASLOG_TMP_ADD(C_);

        if (beta!=BETA(0)) {
            C_ = C;
        }
        mm(transposeA, transposeB, alpha, A, B, beta, C_);
        C = C_;

        FLENS_BLASLOG_TMP_REMOVE(C_, C);
        return;
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_GEMM(transposeA, transposeB, alpha, A, B, beta, C);

#   ifdef HAVE_CXXBLAS_GEMM_STRASSEN
    cxxblas::gemm_strassen(C.order(),
                           transposeA, transposeB,
                           C.numRows(),
                           C.numCols(),
                           k,
                           alpha,
                           A.data(), A.leadingDimension(),
                           B.data(), B.leadingDimension(),
                           beta,
                           C.data(), C.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

} } } // namespace extensions, blas, flens

#endif // PLAYGROUND_FLENS_BLASEXTENSIONS_LEVEL3_MM_TCC
