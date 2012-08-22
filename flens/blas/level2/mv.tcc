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

#include <flens/blas/closures/debugclosure.h>
#include <flens/blas/level2/level2.h>
#include <flens/typedefs.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//== GeneralMatrix - Vector products ===========================================

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

        typename Result<VectorX>::Type  _x;
        FLENS_BLASLOG_TMP_ADD(_x);
        _x = x;

        mv(transpose, alpha, A, _x, beta, y);

        FLENS_BLASLOG_TMP_REMOVE(_x, x);
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

    typedef typename GbMatrix<MA>::IndexType IndexType;
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
        typename Result<DenseVector<VX> >::Type _x;
        FLENS_BLASLOG_TMP_ADD(_x);
        _x = x;

        mv(transpose, alpha, A, _x, beta, y);

        FLENS_BLASLOG_TMP_REMOVE(_x, x);
        return;
    }
#   endif

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
    ASSERT(x.stride()==1);
    ASSERT(y.stride()==1);
#   else
//
//  If x and y are identical an temporary is needed if we want to use mv
//
    if (DEBUGCLOSURE::identical(x, y) || (x.stride()!=1)) {
        typedef typename RemoveRef<VX>::Type  VectorX;

        typename Result<VectorX>::Type  _x;
        FLENS_BLASLOG_TMP_ADD(_x);
        _x = x;

        mv(trans, alpha, A, _x, beta, y);

        FLENS_BLASLOG_TMP_REMOVE(_x, x);
        return;
    }
    if (y.stride()!=1) {
        typedef typename RemoveRef<VY>::Type  VectorY;

        typename Result<VectorY>::Type  _y;
        FLENS_BLASLOG_TMP_ADD(_y);
        _y = y;

        mv(trans, alpha, A, x, beta, _y);

        FLENS_BLASLOG_TMP_REMOVE(_y, y);
        return;
     }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_GEMV(trans, alpha, A, x, beta, y);

#   ifdef HAVE_CXXBLAS_GECRSMV
    cxxblas::gecrsmv(trans,
                     A.numRows(), A.numCols(),
                     alpha,
                     A.engine().values().data() - A.firstRow(),
                     A.engine().rows().data() - A.firstRow(),
                     A.engine().cols().data() - A.firstCol(),
                     x.data() - A.firstCol(),
                     beta,
                     y.data() - A.firstRow());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
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
    typedef typename RemoveRef<MA>::Type   MatrixA;

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

//== SymmetricMatrix - Vector products =========================================

//-- sbmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsSbMatrix<MA>::value
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
        FLENS_BLASLOG_RESIZE_VECTOR(y, A.dim());
        y.resize(A.dim(), 0);
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
        FLENS_BLASLOG_RESIZE_VECTOR(y, A.dim());
        y.resize(A.dim(), 0);
    }

#   ifdef HAVE_CXXBLAS_SYCRSMV
    cxxblas::sycrsmv(A.upLo(),
                     A.dim(),
                     alpha,
                     A.engine().values().data() - A.firstRow(),
                     A.engine().rows().data() - A.firstRow(),
                     A.engine().cols().data() - A.firstCol(),
                     x.data() - A.firstCol(),
                     beta,
                     y.data() - A.firstRow());
#   else
    ASSERT(0);
#   endif
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
    ASSERT(ADDRESS(y)!=ADDRESS(x));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==BETA(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        y.resize(A.dim(), 0);
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

//-- hemv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsHeMatrix<MA>::value
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

#   ifdef HAVE_CXXBLAS_HEMV
    typedef typename RemoveRef<MA>::Type   MatrixA;

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
mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y)
{
    ASSERT(ADDRESS(y)!=ADDRESS(x));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==BETA(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        y.resize(A.dim(), 0);
    }

#   ifdef HAVE_CXXBLAS_HPMV
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

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL2_MV_TCC