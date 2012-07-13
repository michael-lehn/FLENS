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

#ifndef FLENS_BLAS_LEVEL3_MM_TCC
#define FLENS_BLAS_LEVEL3_MM_TCC

#include <flens/typedefs.h>

namespace flens { namespace blas {

//== product type: GeneralMatrix - GeneralMatrix products

//-- gemm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(Transpose transposeA, Transpose transposeB,
   const ALPHA &alpha,
   const GeMatrix<MA> &A, const GeMatrix<MB> &B,
   const BETA &beta,
   GeMatrix<MC> &C)
{
    const bool noTransA = (transposeA==NoTrans || transposeA==Conj);
    const bool noTransB = (transposeB==NoTrans || transposeB==Conj);

#   ifndef NDEBUG
    int kA = (noTransA) ? A.numCols() : A.numRows();
    int kB = (noTransB) ? B.numRows() : B.numCols();
    ASSERT(kA==kB);
#   endif

    typedef typename GeMatrix<MC>::IndexType IndexType;
    IndexType m = (noTransA) ? A.numRows() : A.numCols();
    IndexType n = (noTransB) ? B.numCols() : B.numRows();
    IndexType k = (noTransA) ? A.numCols() : A.numRows();

    if (MC::order!=MA::order) {
        transposeA = Transpose(transposeA ^ Trans);
    }
    if (MC::order!=MB::order) {
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
        typename GeMatrix<MC>::NoView _C;
        FLENS_BLASLOG_TMP_ADD(_C);

        if (beta!=BETA(0)) {
            _C = C;
        }
        mm(transposeA, transposeB, alpha, A, B, beta, _C);
        C = _C;

        FLENS_BLASLOG_TMP_REMOVE(_C, C);
        return;
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_GEMM(transposeA, transposeB, alpha, A, B, beta, C);

#   ifdef HAVE_CXXBLAS_GEMM
    cxxblas::gemm(MC::order,
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

//== product type: TriangularMatrix - GeneralMatrix products

//-- trmm
template <typename ALPHA, typename MA, typename MB>
void
mm(Side side,
   Transpose transA, const ALPHA &alpha, const TrMatrix<MA> &A,
   GeMatrix<MB> &B)
{
#   ifndef NDEBUG
    ASSERT(MB::order==MA::order);
    if (side==Left) {
        assert(A.dim()==B.numRows());
    } else {
        assert(B.numCols()==A.dim());
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_TRMM(side, transA, alpha, A, B);

#   ifdef HAVE_CXXBLAS_TRMM
    cxxblas::trmm(MB::order, side,
                  A.upLo(), transA, A.diag(),
                  B.numRows(), B.numCols(),
                  alpha,
                  A.data(), A.leadingDimension(),
                  B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}


//== product type: SymmetricMatrix - GeneralMatrix products

//-- symm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(Side side,
   const ALPHA &alpha, const SyMatrix<MA> &A, const GeMatrix<MB> &B,
   const BETA &beta, GeMatrix<MC> &C)
{
#   ifndef NDEBUG
    ASSERT(MC::order==MB::order);
    if (side==Left) {
        ASSERT(A.dim()==B.numRows());
    } else {
        ASSERT(B.numCols()==A.dim());
    }
#   endif

    StorageUpLo upLo = (MC::order==MA::order)
                     ? A.upLo()
                     : StorageUpLo(! A.upLo());

    typedef typename GeMatrix<MC>::IndexType IndexType;
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
        typename GeMatrix<MC>::NoView _C;
        FLENS_BLASLOG_TMP_ADD(_C);

        if (beta!=BETA(0)) {
            _C = C;
        }
        mm(side, alpha, A, B, beta, _C);
        C = _C;

        FLENS_BLASLOG_TMP_REMOVE(_C, C);
        return;
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SYMM(side, alpha, A, B, beta, C);

#   ifdef HAVE_CXXBLAS_SYMM
    cxxblas::symm(MC::order, side,
                  upLo,
                  C.numRows(), C.numCols(),
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

//== product type: HermitianMatrix - GeneralMatrix products

//-- hemm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(Side side,
   const ALPHA &alpha, const HeMatrix<MA> &A, const GeMatrix<MB> &B,
   const BETA &beta, GeMatrix<MC> &C)
{
#   ifndef NDEBUG
    ASSERT(MC::order==MB::Oorder);
    if (side==Left) {
        ASSERT(A.dim()==B.numRows());
    } else {
        ASSERT(B.numCols()==A.dim());
    }
#   endif

    StorageUpLo upLo = (MC::order==MA::order)
                     ? A.upLo()
                     : StorageUpLo(! A.upLo());

    typedef typename GeMatrix<MC>::IndexType IndexType;
    IndexType m = (side==Left) ? A.dim() : B.numRows();
    IndexType n = (side==Left) ? B.numCols() : A.dim();
 
    ASSERT((beta==static_cast<BETA>(0)) || (C.numRows()==m));
    ASSERT((beta==static_cast<BETA>(0)) || (C.numCols()==n));
 
    if ((C.numRows()!=m) || (C.numCols()!=n)) {
        C.resize(m,n);
    }

#   ifdef HAVE_CXXBLAS_HEMM
    cxxblas::hemm(MC::order, side,
                  upLo,
                  C.numRows(), C.numCols(),
                  alpha,
                  A.data(), A.leadingDimension(),
                  B.data(), B.leadingDimension(),
                  beta,
                  C.data(), C.leadingDimension());
#   else
    ASSERT(0);
#   endif
}

//== product type: GeneralBandedMatrix - GeneralMatrix products

//-- gbmm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(Transpose transposeA, Transpose transposeB,
   const ALPHA &alpha,
   const GbMatrix<MA> &A, const GeMatrix<MB> &B,
   const BETA &beta,
   GeMatrix<MC> &C)
{
    const bool noTransA = (transposeA==NoTrans || transposeA==Conj);
    const bool noTransB = (transposeB==NoTrans || transposeB==Conj);

#   ifndef NDEBUG
    int kA = (noTransA) ? A.numCols() : A.numRows();
    int kB = (noTransB) ? B.numRows() : B.numCols();
    ASSERT(kA==kB);
#   endif

    typedef typename GeMatrix<MC>::IndexType IndexType;
    IndexType m = (noTransA) ? A.numRows() : A.numCols();
    IndexType n = (noTransB) ? B.numCols() : B.numRows();
    IndexType l = (noTransA) ? A.numCols() : A.numRows();

    if (MC::order!=MA::order) {
        transposeA = Transpose(transposeA ^ Trans);
    }
    if (MC::order!=MB::order) {
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
        C.resize(m, n);
    }
    
#   ifndef FLENS_DEBUG_CLOSURES
    ASSERT(!DEBUGCLOSURE::identical(B, C));
#   else
//
//  If B is identical with C we copy C into a temporary first.  Then
//  we compute the matrix-matrix product and afterwards copy the result into C.
//
    if (DEBUGCLOSURE::identical(B, C)) {
        typename GeMatrix<MC>::NoView _C;
        FLENS_BLASLOG_TMP_ADD(_C);

        if (beta!=BETA(0)) {
            _C = C;
        }
        mm(transposeA, transposeB, alpha, A, B, beta, _C);
        C = _C;

        FLENS_BLASLOG_TMP_REMOVE(_C, C);
        return;
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_GBMM(transposeA, transposeB, alpha, A, B, beta, C);

#   ifdef HAVE_CXXBLAS_GBMM
    cxxblas::gbmm(MC::order, cxxblas::Side::Left, 
                  transposeA, transposeB,
                  A.numRows(), A.numCols(),
                  A.numSubDiags(), A.numSuperDiags(),
                  C.numCols(),
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

//== product type: GeneralMatrix - GeneralBandedMatrix products

//-- gemm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(Transpose transposeA, Transpose transposeB,
   const ALPHA &alpha,
   const GeMatrix<MA> &A, const GbMatrix<MB> &B,
   const BETA &beta,
   GeMatrix<MC> &C)
{
    const bool noTransA = (transposeA==NoTrans || transposeA==Conj);
    const bool noTransB = (transposeB==NoTrans || transposeB==Conj);

#   ifndef NDEBUG
    int kA = (noTransA) ? A.numCols() : A.numRows();
    int kB = (noTransB) ? B.numRows() : B.numCols();
    ASSERT(kA==kB);
#   endif

    typedef typename GeMatrix<MC>::IndexType IndexType;
    IndexType m = (noTransA) ? A.numRows() : A.numCols();
    IndexType n = (noTransB) ? B.numCols() : B.numRows();
    IndexType l = (noTransA) ? A.numCols() : A.numRows();

    if (MC::order!=MA::order) {
        transposeA = Transpose(transposeA ^ Trans);
    }
    if (MC::order!=MB::order) {
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
        C.resize(m, n);
    }
    
#   ifndef FLENS_DEBUG_CLOSURES
    ASSERT(!DEBUGCLOSURE::identical(A, C));
#   else
//
//  If B is identical with C we copy C into a temporary first.  Then
//  we compute the matrix-matrix product and afterwards copy the result into C.
//
    if (DEBUGCLOSURE::identical(A, C)) {
        typename GeMatrix<MC>::NoView _C;
        FLENS_BLASLOG_TMP_ADD(_C);

        if (beta!=BETA(0)) {
            _C = C;
        }
        mm(transposeA, transposeB, alpha, A, B, beta, _C);
        C = _C;

        FLENS_BLASLOG_TMP_REMOVE(_C, C);
        return;
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_GBMM(transposeA, transposeB, alpha, A, B, beta, C);

#   ifdef HAVE_CXXBLAS_GBMM
    cxxblas::gbmm(MC::order, cxxblas::Side::Right, 
                  transposeB, transposeA,
                  A.numRows(), A.numCols(),
                  B.numSubDiags(), B.numSuperDiags(),
                  C.numCols(),
                  alpha,
                  B.data(), B.leadingDimension(),
                  A.data(), A.leadingDimension(),
                  beta,
                  C.data(), C.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}


//== product type: HermitianMatrix - GeneralMatrix products

//-- hbmm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(Side side,
   const ALPHA &alpha, const HbMatrix<MA> &A, const GeMatrix<MB> &B,
   const BETA &beta, GeMatrix<MC> &C)
{
#   ifndef NDEBUG
    ASSERT(MC::order==MB::order);
    if (side==Left) {
        ASSERT(A.dim()==B.numRows());
    } else {
        ASSERT(B.numCols()==A.dim());
    }
#   endif

    StorageUpLo upLo = (MC::order==MA::order)
                     ? A.upLo()
                     : StorageUpLo(! A.upLo());

    typedef typename GeMatrix<MC>::IndexType IndexType;
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
        C.resize(m, n);
    }

#   ifndef FLENS_DEBUG_CLOSURES
    ASSERT(!DEBUGCLOSURE::identical(B, C));
#   else
//
//  If A or B is identical with C we copy C into a temporary first.  Then
//  we compute the matrix-matrix product and afterwards copy the result into C.
//
    if (DEBUGCLOSURE::identical(B, C)) {
        typename GeMatrix<MC>::NoView _C;
        FLENS_BLASLOG_TMP_ADD(_C);

        if (beta!=BETA(0)) {
            _C = C;
        }
        mm(side, alpha, A, B, beta, _C);
        C = _C;

        FLENS_BLASLOG_TMP_REMOVE(_C, C);
        return;
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SYMM(side, alpha, A, B, beta, C);

#   ifdef HAVE_CXXBLAS_HBMM
    cxxblas::hbmm(MC::order, side,
                  upLo,
                  C.numRows(), A.numOffDiags(), C.numCols(),
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
//== product type: SymmetricBandedMatrix - GeneralMatrix products

//-- sbmm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(Side side,
   const ALPHA &alpha, const SbMatrix<MA> &A, const GeMatrix<MB> &B,
   const BETA &beta, GeMatrix<MC> &C)
{
#   ifndef NDEBUG
    ASSERT(MC::order==MB::order);
    if (side==Left) {
        ASSERT(A.dim()==B.numRows());
    } else {
        ASSERT(B.numCols()==A.dim());
    }
#   endif

    StorageUpLo upLo = (MC::order==MA::order)
                     ? A.upLo()
                     : StorageUpLo(! A.upLo());

    typedef typename GeMatrix<MC>::IndexType IndexType;
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
        C.resize(m, n);
    }

#   ifndef FLENS_DEBUG_CLOSURES
    ASSERT(!DEBUGCLOSURE::identical(B, C));
#   else
//
//  If A or B is identical with C we copy C into a temporary first.  Then
//  we compute the matrix-matrix product and afterwards copy the result into C.
//
    if (DEBUGCLOSURE::identical(B, C)) {
        typename GeMatrix<MC>::NoView _C;
        FLENS_BLASLOG_TMP_ADD(_C);

        if (beta!=BETA(0)) {
            _C = C;
        }
        mm(side, alpha, A, B, beta, _C);
        C = _C;

        FLENS_BLASLOG_TMP_REMOVE(_C, C);
        return;
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SYMM(side, alpha, A, B, beta, C);

#   ifdef HAVE_CXXBLAS_SBMM
    cxxblas::sbmm(MC::order, side,
                  upLo,
                  C.numRows(), A.numOffDiags(), C.numCols(),
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


//== product type: TriangularBandedMatrix - GeneralMatrix products

//-- tbmm
template <typename ALPHA, typename MA, typename MB>
void
mm(Side side,
   Transpose transA, const ALPHA &alpha, const TbMatrix<MA> &A,
   GeMatrix<MB> &B)
{
#   ifndef NDEBUG
    ASSERT(MB::order==MA::order);
    if (side==Left) {
        assert(A.dim()==B.numRows());
    } else {
        assert(B.numCols()==A.dim());
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_TRMM(side, transA, alpha, A, B);

#   ifdef HAVE_CXXBLAS_TBMM
    cxxblas::tbmm(MB::order, side,
                  A.upLo(), transA, A.diag(),
                  B.numRows(), B.numCols(),
                  A.numOffDiags(),
                  alpha,
                  A.data(), A.leadingDimension(),
                  B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}


//== Forwarding ================================================================

//-- GeneralMatrix - GeneralMatrix products
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
typename RestrictTo<IsGeneralMatrix<MA>::value &&
                    IsGeneralMatrix<MB>::value &&
                   !IsClosure<MA>::value &&
                   !IsClosure<MB>::value &&
                    IsSame<MC, typename MC::Impl>::value,
         void>::Type
mm(Transpose transA, Transpose transB, const ALPHA &alpha,
   const MA &A, const MB &B, const BETA &beta, MC &&C)
{
    CHECKPOINT_ENTER;
    mm(transA, transB, alpha, A, B, beta, C);
    CHECKPOINT_LEAVE;
}

//-- TriangularMatrix - GeneralMatrix products
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsTriangularMatrix<MA>::value &&
                    IsGeneralMatrix<MB>::value &&
                   !IsClosure<MA>::value &&
                    IsSame<MB, typename MB::Impl>::value,
         void>::Type
mm(Side side, Transpose transA, const ALPHA &alpha, const MA &A, MB &&B)
{
    CHECKPOINT_ENTER;
    mm(side, transA, alpha, A, B);
    CHECKPOINT_LEAVE;
}


//-- SymmetricMatrix - GeneralMatrix products
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
typename RestrictTo<IsSymmetricMatrix<MA>::value &&
                    IsGeneralMatrix<MB>::value &&
                   !IsClosure<MA>::value &&
                   !IsClosure<MB>::value &&
                    IsSame<MC, typename MC::Impl>::value,
         void>::Type
mm(Side side, const ALPHA &alpha, const MA &A, const MB &B,
   const BETA &beta, MC &&C)
{
    CHECKPOINT_ENTER;
    mm(side, alpha, A, B, beta, C);
    CHECKPOINT_LEAVE;
}

//-- HermitianMatrix - GeneralMatrix products
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
typename RestrictTo<IsHermitianMatrix<MA>::value &&
                    IsGeneralMatrix<MB>::value &&
                   !IsClosure<MA>::value &&
                   !IsClosure<MB>::value &&
                    IsSame<MC, typename MC::Impl>::value,
         void>::Type
mm(Side side, const ALPHA &alpha, const MA &A, const MB &B,
   const BETA &beta, MC &&C)
{
    CHECKPOINT_ENTER;
    mm(side, alpha, A, B, beta, C);
    CHECKPOINT_LEAVE;
}


} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL3_MM_TCC
