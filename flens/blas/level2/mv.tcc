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

#ifndef FLENS_BLAS_LEVEL2_MV_TCC
#define FLENS_BLAS_LEVEL2_MV_TCC 1

#include <flens/blas/closures/closures.h>
#include <flens/blas/level1/level1.h>
#include <flens/blas/level2/level2.h>
#include <flens/typedefs.h>
#include <ulmblas/cxxblas.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//== GeneralMatrix - Vector products ===========================================

//--gbmv
template <typename Alpha, typename MA, typename VX, typename Beta, typename VY>
typename RestrictTo<IsGbMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(Transpose transposeA, const Alpha &alpha, const MA &A, const VX &x,
   const Beta &beta, VY &&y)
{
    const bool colMajorA = (A.order()==ColMajor);
    const bool transA    = (transposeA==Trans || transposeA==ConjTrans);
    const bool conjA     = (transposeA==Conj || transposeA==ConjTrans);


#   ifndef NDEBUG
    if (!transA) {
        ASSERT(x.length()==A.numCols());
    } else {
        ASSERT(x.length()==A.numRows());
    }
#   endif

    typedef typename RemoveRef<MA>::Type  MatrixA;
    typedef typename MatrixA::IndexType   IndexType;

    IndexType yLength = (!transA) ? A.numRows() : A.numCols();

    ASSERT(beta==Beta(0) || y.length()==yLength || y.length()==0);

    if (y.length()!=yLength) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        FLENS_BLASLOG_RESIZE_VECTOR(y, yLength);
        y.resize(yLength, y.firstIndex(), Zero);
    }

    ASSERT(!DEBUGCLOSURE::identical(x, y));

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_GBMV(transposeA, alpha, A, x, beta, y);

    cxxblas::gbmv(A.numRows(), A.numCols(), A.numSubDiags(), A.numSuperDiags(),
                  alpha,
                  colMajorA, transA, conjA, A.data(), A.leadingDimension(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- gemv
template <typename Alpha, typename MA, typename VX, typename Beta, typename VY>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(Transpose transpose, const Alpha &alpha, const MA &A, const VX &x,
   const Beta &beta, VY &&y)
{
    const bool transA = (transpose==Trans || transpose==ConjTrans);
    const bool conjA  = (transpose==Conj || transpose==ConjTrans);

#   ifndef NDEBUG
    if (!transA) {
        ASSERT(x.length()==A.numCols());
    } else {
        ASSERT(x.length()==A.numRows());
    }
#   endif

    typedef typename RemoveRef<MA>::Type  MatrixA;
    typedef typename MatrixA::IndexType   IndexType;

    const IndexType yLength = (!transA) ? A.numRows() : A.numCols();

    ASSERT(beta==Beta(0) || y.length()==yLength || y.length()==0);

    if (y.length()!=yLength) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        FLENS_BLASLOG_RESIZE_VECTOR(y, yLength);
        y.resize(yLength, y.firstIndex(), Zero);
    }


#   ifndef FLENS_DEBUG_CLOSURES
    ASSERT(!DEBUGCLOSURE::identical(x, y));
#   else
//
//  If x and y are identical an temporary is needed if we want to use mv
//
    if (DEBUGCLOSURE::identical(x, y)) {
        typedef typename RemoveRef<VX>::Type  VectorX;

        typename Result<VectorX>::Type  x_;
        FLENS_BLASLOG_TMP_ADD(x_);
        x_ = x;

        mv(transpose, alpha, A, x_, beta, y);

        FLENS_BLASLOG_TMP_REMOVE(x_, x);
        return;
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_GEMV(transpose, alpha, A, x, beta, y);

    cxxblas::gemv(A.numRows(), A.numCols(),
                  alpha,
                  transA, conjA,
                  A.data(), A.strideRow(), A.strideCol(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//== HermitianMatrix - Vector products =========================================

//-- hbmv
template <typename Alpha, typename MA, typename VX, typename Beta, typename VY>
typename RestrictTo<IsHbMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const Alpha &alpha, const MA &A, const VX &x, const Beta &beta, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==Beta(0)) || (y.length()==A.dim()));

    const bool colMajorA = (A.order()==ColMajor);
    const bool lowerA    = (A.upLo()==Lower);

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        FLENS_BLASLOG_RESIZE_VECTOR(y, A.dim());
        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

    cxxblas::hbmv(A.dim(), A.numOffDiags(),
                  alpha,
                  colMajorA, lowerA,
                  A.data(), A.leadingDimension(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());
}

//-- hemv
template <typename Alpha, typename MA, typename VX, typename Beta, typename VY>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const Alpha &alpha, const MA &A, const VX &x, const Beta &beta, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==Beta(0)) || (y.length()==A.dim()));

    const bool lowerA    = (A.upLo()==Lower);

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

    cxxblas::hemv(A.dim(),
                  alpha,
                  lowerA, A.data(), A.strideRow(), A.strideCol(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());
}

//-- hpmv
template <typename Alpha, typename MA, typename VX, typename Beta, typename VY>
typename RestrictTo<IsHpMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const Alpha &alpha, const MA &A, const VX &x,
   const Beta &beta, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==Beta(0)) || (y.length()==A.dim()));

    const bool colMajorA = (A.order()==ColMajor);
    const bool lowerA    = (A.upLo()==Lower);

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

    cxxblas::hpmv(A.dim(),
                  alpha,
                  colMajorA, lowerA, A.data(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());
}


//== SymmetricMatrix - Vector products =========================================

//-- sbmv
template <typename Alpha, typename MA, typename VX, typename Beta, typename VY>
typename RestrictTo<IsSbMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const Alpha &alpha, const MA &A, const VX &x, const Beta &beta, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==Beta(0)) || (y.length()==A.dim()));

    const bool colMajorA = (A.order()==ColMajor);
    const bool lowerA    = (A.upLo()==Lower);

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

    cxxblas::sbmv(A.dim(), A.numOffDiags(),
                  alpha,
                  colMajorA, lowerA, A.data(), A.leadingDimension(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());
}

//-- spmv
template <typename Alpha, typename MA, typename VX, typename Beta, typename VY>
typename RestrictTo<IsSpMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const Alpha &alpha, const MA &A, const VX &x, const Beta &beta, VY &&y)
{
    ASSERT(ADDRESS(y)!=ADDRESS(x));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==Beta(0)) || (y.length()==A.dim()));

    const bool colMajorA = (A.order()==ColMajor);
    const bool lowerA    = (A.upLo()==Lower);

    if (y.length()!=A.dim()) {
        y.resize(A.dim(), 0);
    }

    cxxblas::spmv(A.dim(),
                  alpha,
                  colMajorA, lowerA, A.data(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());
}

//-- symv
template <typename Alpha, typename MA, typename VX, typename Beta, typename VY>
typename RestrictTo<IsSyMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const Alpha &alpha, const MA &A, const VX &x, const Beta &beta, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==Beta(0)) || (y.length()==A.dim()));

    const bool lowerA    = (A.upLo()==Lower);

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        FLENS_BLASLOG_RESIZE_VECTOR(y, A.dim());
        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

    cxxblas::symv(A.dim(),
                  alpha,
                  lowerA, A.data(), A.strideRow(), A.strideCol(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());
}

//== TriangularMatrix - Vector products ========================================

//-- tbmv
template <typename MA, typename VX>
typename RestrictTo<IsTbMatrix<MA>::value
                 && IsDenseVector<VX>::value,
         void>::Type
mv(Transpose transposeA, const MA &A, VX &&x)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_TBMV(transposeA, A, x);

    const bool colMajorA = (A.order()==ColMajor);
    const bool lowerA    = (A.upLo()==Lower);
    const bool transA    = (transposeA==Trans || transposeA==ConjTrans);
    const bool conjA     = (transposeA==Conj || transposeA==ConjTrans);
    const bool unitDiagA = (A.diag()==Unit);

#   ifndef NDEBUG
    if (NoTrans) {
        ASSERT(x.length()==A.numCols());
    } else {
        ASSERT(x.length()==A.numRows());
    }
#   endif

    cxxblas::tbmv(A.dim(), A.numOffDiags(),
                  colMajorA, lowerA, transA, conjA, unitDiagA,
                  A.data(), A.leadingDimension(),
                  x.data(), x.stride());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- trmv
template <typename MA, typename VX>
typename RestrictTo<IsTrMatrix<MA>::value
                 && IsDenseVector<VX>::value,
         void>::Type
mv(Transpose transposeA, const MA &A, VX &&x)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_TRMV(transposeA, A, x);

    ASSERT(x.length()==A.dim());

    const bool lowerA    = (A.upLo()==Lower);
    const bool transA    = (transposeA==Trans || transposeA==ConjTrans);
    const bool conjA     = (transposeA==Conj || transposeA==ConjTrans);
    const bool unitDiagA = (A.diag()==Unit);

    cxxblas::trmv(A.dim(),
                  lowerA, transA, conjA, unitDiagA,
                  A.data(), A.strideRow(), A.strideCol(),
                  x.data(), x.stride());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- tpmv
template <typename MA, typename VX>
typename RestrictTo<IsTpMatrix<MA>::value
                 && IsDenseVector<VX>::value,
         void>::Type
mv(Transpose transposeA, const MA &A, VX &&x)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_TPMV(transposeA, A, x);

    ASSERT(x.length()==A.dim());

    const bool colMajorA = (A.order()==ColMajor);
    const bool lowerA    = (A.upLo()==Lower);
    const bool transA    = (transposeA==Trans || transposeA==ConjTrans);
    const bool conjA     = (transposeA==Conj || transposeA==ConjTrans);
    const bool unitDiagA = (A.diag()==Unit);

    cxxblas::tpmv(A.dim(),
                  colMajorA, lowerA, transA, conjA, unitDiagA,
                  A.data(),
                  x.data(), x.stride());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL2_MV_TCC
