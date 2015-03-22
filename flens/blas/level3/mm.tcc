/*
 *   Copyright (c) 2009,2015 Michael Lehn
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

#ifndef FLENS_BLAS_LEVEL3_MM_TCC
#define FLENS_BLAS_LEVEL3_MM_TCC

#include <flens/blas/closures/closures.h>
#include <flens/blas/level3/level3.h>
#include <flens/typedefs.h>
#include <ulmblas/cxxblas.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//== GeneralMatrix - GeneralMatrix products ====================================

//-- gemm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsGeMatrix<MB>::value
                 && IsGeMatrix<MC>::value,
         void>::Type
mm(Transpose        transposeA,
   Transpose        transposeB,
   const ALPHA      &alpha,
   const MA         &A,
   const MB         &B,
   const BETA       &beta,
   MC               &&C)
{
    typedef typename RemoveRef<MC>::Type MatrixC;
    typedef typename MatrixC::IndexType  IndexType;

    const bool transA = (transposeA==Trans || transposeA==ConjTrans);
    const bool transB = (transposeB==Trans || transposeB==ConjTrans);

    const bool conjA = (transposeA==Conj || transposeA==ConjTrans);
    const bool conjB = (transposeB==Conj || transposeB==ConjTrans);


#   ifndef NDEBUG
    int kA = (! transA) ? A.numCols() : A.numRows();
    int kB = (! transB) ? B.numRows() : B.numCols();
    ASSERT(kA==kB);
#   endif

    IndexType m = (! transA) ? A.numRows() : A.numCols();
    IndexType n = (! transB) ? B.numCols() : B.numRows();
    IndexType k = (! transA) ? A.numCols() : A.numRows();

#   ifndef NDEBUG
    if (beta!=BETA(0)) {
        if (C.numRows()!=0 && C.numCols()!=0) {
            if (C.numRows()!=m || C.numCols()!=n) {
                std::cerr << "A = " << A << std::endl;
                std::cerr << "transposeA = " << transposeA << std::endl;
                std::cerr << "B = " << B << std::endl;
                std::cerr << "transposeB = " << transposeB << std::endl;
                std::cerr << "C = " << C << std::endl;
            }
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

    cxxblas::gemm(C.numRows(), C.numCols(), k,
                  alpha,
                  transA, conjA, A.data(), A.strideRow(), A.strideCol(),
                  transB, conjB, B.data(), B.strideRow(), B.strideCol(),
                  beta,
                  C.data(), C.strideRow(), C.strideCol());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//== HermitianMatrix - GeneralMatrix products ==================================

//-- hemm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsGeMatrix<MB>::value
                 && IsGeMatrix<MC>::value,
         void>::Type
mm(Side             side,
   const ALPHA      &alpha,
   const MA         &A,
   const MB         &B,
   const BETA       &beta,
   MC               &&C)
{
    typedef typename RemoveRef<MC>::Type MatrixC;
    typedef typename MatrixC::IndexType  IndexType;

#   ifndef NDEBUG
    if (side==Left) {
        ASSERT(A.dim()==B.numRows());
    } else {
        ASSERT(B.numCols()==A.dim());
    }
#   endif

    IndexType m = (side==Left) ? A.dim() : B.numRows();
    IndexType n = (side==Left) ? B.numCols() : A.dim();

    ASSERT((beta==static_cast<BETA>(0)) || (C.numRows()==m));
    ASSERT((beta==static_cast<BETA>(0)) || (C.numCols()==n));

    if ((C.numRows()!=m) || (C.numCols()!=n)) {
        ASSERT(C.numRows()==0 && C.numCols()==0);
        C.resize(m,n);
    }

    cxxblas::hemm(side==Left,
                  C.numRows(), C.numCols(),
                  alpha,
                  A.upLo()==Lower, A.data(), A.strideRow(), A.strideCol(),
                  B.data(), B.strideRow(), B.strideCol(),
                  beta,
                  C.data(), C.strideRow(), C.strideCol());
}

//== SymmetricMatrix - GeneralMatrix products ==================================

//--symm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
typename RestrictTo<IsSyMatrix<MA>::value
                 && IsGeMatrix<MB>::value
                 && IsGeMatrix<MC>::value,
         void>::Type
mm(Side             side,
   const ALPHA      &alpha,
   const MA         &A,
   const MB         &B,
   const BETA       &beta,
   MC               &&C)
{
    typedef typename RemoveRef<MC>::Type MatrixC;
    typedef typename MatrixC::IndexType  IndexType;

#   ifndef NDEBUG
    if (side==Left) {
        ASSERT(A.dim()==B.numRows());
    } else {
        ASSERT(B.numCols()==A.dim());
    }
#   endif

    IndexType m = (side==Left) ? A.dim() : B.numRows();
    IndexType n = (side==Left) ? B.numCols() : A.dim();

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
        mm(side, alpha, A, B, beta, C_);
        C = C_;

        FLENS_BLASLOG_TMP_REMOVE(C_, C);
        return;
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SYMM(side, alpha, A, B, beta, C);

    cxxblas::symm(side==Left,
                  C.numRows(), C.numCols(),
                  alpha,
                  A.upLo()==Lower, A.data(), A.strideRow(), A.strideCol(),
                  B.data(), B.strideRow(), B.strideCol(),
                  beta,
                  C.data(), C.strideRow(), C.strideCol());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}


//== TriangularMatrix - GeneralMatrix products =================================

//-- trmm
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsTrMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
mm(Side             side,
   Transpose        transposeA,
   const ALPHA      &alpha,
   const MA         &A,
   MB               &&B)
{
#   ifndef NDEBUG
    if (side==Left) {
        assert(A.dim()==B.numRows());
    } else {
        assert(B.numCols()==A.dim());
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_TRMM(side, transA, alpha, A, B);

    const bool lowerA    = (A.upLo()==Lower);
    const bool transA    = (transposeA==Trans || transposeA==ConjTrans);
    const bool conjA     = (transposeA==Conj || transposeA==ConjTrans);
    const bool unitDiagA = (A.diag()==Unit);

    cxxblas::trmm(side==Left,
                  B.numRows(), B.numCols(),
                  alpha,
                  lowerA, transA, conjA, unitDiagA,
                  A.data(), A.strideRow(), A.strideCol(),
                  B.data(), B.strideRow(), B.strideCol());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL3_MM_TCC
