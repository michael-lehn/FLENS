/*
 *   Copyright (c) 2009, Michael Lehn
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

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//== GeneralMatrix - Vector products ===========================================

//--diagmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsDiagMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(Transpose transpose, const ALPHA &alpha, const MA &A, const VX &x,
   const BETA &beta, VY &&y)
{
    typedef typename RemoveRef<MA>::Type MatrixA;
    typedef typename MatrixA::IndexType  IndexType;

    const bool noTrans = (transpose==NoTrans || transpose==Conj);

#   ifndef NDEBUG
    if (noTrans) {
        ASSERT(x.length()==A.numCols());
    } else {
        ASSERT(x.length()==A.numRows());
    }
#   endif

    IndexType yLength = noTrans ? A.numRows()
                                : A.numCols();

    ASSERT(beta==BETA(0) || y.length()==yLength || y.length()==0);

    if (y.length()!=yLength) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        FLENS_BLASLOG_RESIZE_VECTOR(y, yLength);
        y.resize(yLength, y.firstIndex(), Zero);
    }

    if (DEBUGCLOSURE::identical(x, y)) {
        ASSERT( beta==BETA(0) );

        blas::scal(alpha, y);

        FLENS_BLASLOG_SETTAG("--> ");
        FLENS_BLASLOG_BEGIN_TBMV(transpose, A, x);
#       ifdef HAVE_CXXBLAS_TBMV
        cxxblas::tbmv(ColMajor, Upper,
                      transpose, NonUnit,
                      A.dim(), IndexType(0),
                      A.data(), A.leadingDimension(),
                      y.data(), y.stride());
#       else
        ASSERT(0);
#       endif
        FLENS_BLASLOG_END;
        FLENS_BLASLOG_UNSETTAG;
        return;
    }

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_GBMV(transpose, alpha, A, x, beta, y);

#   ifdef HAVE_CXXBLAS_GBMV
    cxxblas::gbmv(ColMajor,
                  transpose,
                  A.numRows(), A.numCols(),
                  IndexType(0), IndexType(0),
                  alpha,
                  A.data(), A.leadingDimension(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//--gbmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsGbMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(Transpose transpose, const ALPHA &alpha, const MA &A, const VX &x,
   const BETA &beta, VY &&y)
{

    const bool noTrans = (transpose==NoTrans || transpose==Conj);

#   ifndef NDEBUG
    if (noTrans) {
        ASSERT(x.length()==A.numCols());
    } else {
        ASSERT(x.length()==A.numRows());
    }
#   endif

    typedef typename RemoveRef<MA>::Type  MatrixA;
    typedef typename MatrixA::IndexType   IndexType;
    IndexType yLength = noTrans ? A.numRows()
                                : A.numCols();

    ASSERT(beta==BETA(0) || y.length()==yLength || y.length()==0);

    if (y.length()!=yLength) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        FLENS_BLASLOG_RESIZE_VECTOR(y, yLength);
        y.resize(yLength, y.firstIndex(), Zero);
    }

    ASSERT(!DEBUGCLOSURE::identical(x, y));

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_GBMV(transpose, alpha, A, x, beta, y);

#   ifdef HAVE_CXXBLAS_GBMV
    cxxblas::gbmv(A.order(),
                  transpose,
                  A.numRows(), A.numCols(),
                  A.numSubDiags(), A.numSuperDiags(),
                  alpha,
                  A.data(), A.leadingDimension(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- geccsmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsGeCCSMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(Transpose trans, const ALPHA &alpha, const MA &A, const VX &x,
   const BETA &beta, VY &&y)
{
    const bool noTrans = (trans==NoTrans || trans==Conj);

#   ifndef NDEBUG
    if (noTrans) {
        ASSERT(x.length()==A.numCols());
    } else {
        ASSERT(x.length()==A.numRows());
    }
#   endif

    typedef typename RemoveRef<MA>::Type  MatrixA;
    typedef typename MatrixA::IndexType   IndexType;
    IndexType yLength = noTrans ? A.numRows()
                                : A.numCols();

    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(beta==BETA(0) || y.length()==yLength || y.length()==0);

    if (y.length()!=yLength) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(yLength, y.firstIndex(), Zero);
    }

//  Sparse BLAS only supports this case:
    ASSERT(x.stride()==1);
    ASSERT(y.stride()==1);

#   ifdef HAVE_CXXBLAS_GECRSMV
    cxxblas::gecrsmv(Transpose(trans^Trans),
                     A.numCols(), A.numRows(),
                     alpha,
                     A.engine().values().data(),
                     A.engine().cols().data(),
                     A.engine().rows().data(),
                     x.data(),
                     beta,
                     y.data());
#   else
    ASSERT(0);
#   endif
}

//-- gecrsmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsGeCRSMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(Transpose trans, const ALPHA &alpha, const MA &A, const VX &x,
   const BETA &beta, VY &&y)
{
    const bool noTrans = (trans==NoTrans || trans==Conj);

#   ifndef NDEBUG
    if (noTrans) {
        ASSERT(x.length()==A.numCols());
    } else {
        ASSERT(x.length()==A.numRows());
    }
#   endif

    typedef typename RemoveRef<MA>::Type  MatrixA;
    typedef typename MatrixA::IndexType   IndexType;
    IndexType yLength = noTrans ? A.numRows()
                                : A.numCols();

    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(beta==BETA(0) || y.length()==yLength || y.length()==0);

    if (y.length()!=yLength) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(yLength, y.firstIndex(), Zero);
    }

//  Sparse BLAS only supports this case:
    ASSERT(x.stride()==1);
    ASSERT(y.stride()==1);

#   ifdef HAVE_CXXBLAS_GECRSMV
    cxxblas::gecrsmv(trans,
                     A.numRows(), A.numCols(),
                     alpha,
                     A.engine().values().data(),
                     A.engine().rows().data(),
                     A.engine().cols().data(),
                     x.data(),
                     beta,
                     y.data());
#   else
    ASSERT(0);
#   endif
}

//-- gemv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(Transpose transpose, const ALPHA &alpha, const MA &A, const VX &x,
   const BETA &beta, VY &&y)
{
    const bool noTrans = (transpose==NoTrans || transpose==Conj);

#   ifndef NDEBUG
    if (noTrans) {
        ASSERT(x.length()==A.numCols());
    } else {
        ASSERT(x.length()==A.numRows());
    }
#   endif

    typedef typename RemoveRef<MA>::Type  MatrixA;
    typedef typename MatrixA::IndexType   IndexType;
    IndexType yLength = noTrans ? A.numRows()
                                : A.numCols();

    ASSERT(beta==BETA(0) || y.length()==yLength || y.length()==0);

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

#   ifdef HAVE_CXXBLAS_GEMV
    cxxblas::gemv(A.order(),
                  transpose,
                  A.numRows(), A.numCols(),
                  alpha,
                  A.data(), A.leadingDimension(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- (tiny) gemv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsGeTinyMatrix<MA>::value
                 && IsTinyVector<VX>::value
                 && IsTinyVector<VY>::value,
         void>::Type
mv(Transpose trans, const ALPHA &alpha, const MA &A, const VX &x,
   const BETA &beta, VY &&y)
{
#   ifndef NDEBUG
    if (trans==NoTrans || trans==ConjTrans) {
        ASSERT(A.numCols()==x.length());
        ASSERT(A.numRows()==y.length());
    } else {
        ASSERT(A.numCols()==y.length());
        ASSERT(A.numRows()==x.length());
    }
#   endif

    typedef typename MA::ElementType       TA;
    typedef typename VX::ElementType       TX;
    typedef typename RemoveRef<VY>::Type   VectorY;
    typedef typename VectorY::ElementType  TY;

    const int m    = MA::Engine::numRows;
    const int n    = MA::Engine::numCols;
    const int ldA  = MA::Engine::leadingDimension;
    const int incX = VX::Engine::stride;
    const int incY = VectorY::Engine::stride;


    cxxblas::gemv<m,n,TA,ldA,TX,incX,TY,incY>(trans,
                                              alpha, A.data(), x.data(),
                                              beta, y.data());
}

//== HermitianMatrix - Vector products =========================================

//-- hbmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsHbMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==BETA(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        FLENS_BLASLOG_RESIZE_VECTOR(y, A.dim());
        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

#   ifdef HAVE_CXXBLAS_HBMV
    cxxblas::hbmv(A.order(), A.upLo(),
                  A.dim(), A.numOffDiags(),
                  alpha,
                  A.data(), A.leadingDimension(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());
#   else
    ASSERT(0);
#   endif
}

//-- heccsmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsHeCCSMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y)
{
    ASSERT(ADDRESS(y)!=ADDRESS(x));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==BETA(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

#   ifdef HAVE_CXXBLAS_HECRSMV
    cxxblas::heccsmv(A.upLo(),
                     A.dim(),
                     alpha,
                     A.engine().values().data(),
                     A.engine().rows().data(),
                     A.engine().cols().data(),
                     x.data(),
                     beta,
                     y.data());
#   else
    ASSERT(0);
#   endif
}

//-- hecrsmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsHeCRSMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==BETA(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

#   ifdef HAVE_CXXBLAS_HECRSMV
    cxxblas::hecrsmv(A.upLo(),
                     A.dim(),
                     alpha,
                     A.engine().values().data(),
                     A.engine().rows().data(),
                     A.engine().cols().data(),
                     x.data(),
                     beta,
                     y.data());
#   else
    ASSERT(0);
#   endif
}

//-- hemv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==BETA(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

#   ifdef HAVE_CXXBLAS_HEMV
    cxxblas::hemv(A.order(), A.upLo(),
                  A.dim(),
                  alpha,
                  A.data(), A.leadingDimension(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());
#   else
    ASSERT(0);
#   endif
}

//-- hpmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsHpMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const ALPHA &alpha, const MA &A, const VX &x,
   const BETA &beta, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==BETA(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

#   ifdef HAVE_CXXBLAS_HEMV
    cxxblas::hpmv(A.order(), A.upLo(),
                  A.dim(),
                  alpha,
                  A.data(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());
#   else
    ASSERT(0);
#   endif
}


//== SymmetricMatrix - Vector products =========================================

//-- sbmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsSbMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==BETA(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

#   ifdef HAVE_CXXBLAS_SBMV
    cxxblas::sbmv(A.order(), A.upLo(),
                  A.dim(), A.numOffDiags(),
                  alpha,
                  A.data(), A.leadingDimension(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());
#   else
    ASSERT(0);
#   endif
}

//-- spmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsSpMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y)
{
    ASSERT(ADDRESS(y)!=ADDRESS(x));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==BETA(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        y.resize(A.dim(), 0);
    }

#   ifdef HAVE_CXXBLAS_SPMV
    cxxblas::spmv(A.order(), A.upLo(),
                  A.dim(),
                  alpha,
                  A.data(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());
#   else
    ASSERT(0);
#   endif
}

//-- syccsmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsSyCCSMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==BETA(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

#   ifdef HAVE_CXXBLAS_SYCCSMV
    cxxblas::syccsmv(A.upLo(),
                     A.dim(),
                     alpha,
                     A.engine().values().data(),
                     A.engine().rows().data(),
                     A.engine().cols().data(),
                     x.data(),
                     beta,
                     y.data());
#   else
    ASSERT(0);
#   endif
}

//-- sycrsmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsSyCRSMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==BETA(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

#   ifdef HAVE_CXXBLAS_SYCRSMV
    cxxblas::sycrsmv(A.upLo(),
                     A.dim(),
                     alpha,
                     A.engine().values().data(),
                     A.engine().rows().data(),
                     A.engine().cols().data(),
                     x.data(),
                     beta,
                     y.data());
#   else
    ASSERT(0);
#   endif
}

//-- symv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsSyMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==BETA(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        FLENS_BLASLOG_RESIZE_VECTOR(y, A.dim());
        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

#   ifdef HAVE_CXXBLAS_SYMV
    cxxblas::symv(A.order(), A.upLo(),
                  A.dim(),
                  alpha,
                  A.data(), A.leadingDimension(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());
#   else
    ASSERT(0);
#   endif
}

//== TriangularMatrix - Vector products ========================================

//-- tbmv
template <typename MA, typename VX>
typename RestrictTo<IsTbMatrix<MA>::value
                 && IsDenseVector<VX>::value,
         void>::Type
mv(Transpose trans, const MA &A, VX &&x)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_TBMV(trans, A, x);

#   ifndef NDEBUG
    if (NoTrans) {
        ASSERT(x.length()==A.numCols());
    } else {
        ASSERT(x.length()==A.numRows());
    }
#   endif

#   ifdef HAVE_CXXBLAS_TBMV
    cxxblas::tbmv(A.order(), A.upLo(),
                  trans, A.diag(),
                  A.dim(), A.numOffDiags(),
                  A.data(), A.leadingDimension(),
                  x.data(), x.stride());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- tpmv
template <typename MA, typename VX>
typename RestrictTo<IsTpMatrix<MA>::value
                 && IsDenseVector<VX>::value,
         void>::Type
mv(Transpose trans, const MA &A, VX &&x)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_TPMV(trans, A, x);

    ASSERT(x.length()==A.dim());

#   ifdef HAVE_CXXBLAS_TPMV
    cxxblas::tpmv(A.order(), A.upLo(),
                  trans, A.diag(),
                  A.dim(),
                  A.data(),
                  x.data(), x.stride());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- trmv
template <typename MA, typename VX>
typename RestrictTo<IsTrMatrix<MA>::value
                 && IsDenseVector<VX>::value,
         void>::Type
mv(Transpose trans, const MA &A, VX &&x)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_TRMV(trans, A, x);

    ASSERT(x.length()==A.dim());

#   ifdef HAVE_CXXBLAS_TRMV
    cxxblas::trmv(A.order(), A.upLo(),
                  trans, A.diag(),
                  A.dim(),
                  A.data(), A.leadingDimension(),
                  x.data(), x.stride());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL2_MV_TCC
