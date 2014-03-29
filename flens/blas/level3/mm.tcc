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

#include <flens/blas/closures/closures.h>
#include <flens/blas/level3/level3.h>
#include <flens/typedefs.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//== GeneralMatrix - GeneralMatrix products ====================================


//-- diagmm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
typename RestrictTo<IsDiagMatrix<MA>::value
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
    const bool noTransA = (transposeA==NoTrans || transposeA==Conj);
    const bool noTransB = (transposeB==NoTrans || transposeB==Conj);

    typedef typename RemoveRef<MC>::Type MatrixC;
    typedef typename MatrixC::IndexType  IndexType;

#   ifndef NDEBUG
    int kA = (noTransA) ? A.numCols() : A.numRows();
    int kB = (noTransB) ? B.numRows() : B.numCols();
    ASSERT(kA==kB);
#   endif

    IndexType m = (noTransA) ? A.numRows() : A.numCols();
    IndexType n = (noTransB) ? B.numCols() : B.numRows();

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
        C.resize(m, n);
    }

    if(DEBUGCLOSURE::identical(B, C)) {
        ASSERT( beta==BETA(0) );
        ASSERT( transposeB==NoTrans );

        FLENS_BLASLOG_SETTAG("--> ");
        FLENS_BLASLOG_BEGIN_TBMM(side, transA, alpha, A, B);

#       ifdef HAVE_CXXBLAS_TBMM
        cxxblas::tbmm(C.order(), Left,
                      Upper, transposeA, NonUnit,
                      C.numRows(), C.numCols(),
                      IndexType(0),
                      alpha,
                      A.data(), A.leadingDimension(),
                      C.data(), C.leadingDimension());
#       else
        ASSERT(0);
#       endif

        FLENS_BLASLOG_END;
        FLENS_BLASLOG_UNSETTAG;

        return;
    }

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_GBMM(transposeA, transposeB, alpha, A, B, beta, C);

#   ifdef HAVE_CXXBLAS_GBMM
    cxxblas::gbmm(C.order(), cxxblas::Side::Left,
                  transposeA, transposeB,
                  A.numRows(), A.numCols(),
                  IndexType(0), IndexType(0),
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

//-- diagmm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsDiagMatrix<MB>::value
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
    const bool noTransA = (transposeA==NoTrans || transposeA==Conj);
    const bool noTransB = (transposeB==NoTrans || transposeB==Conj);

    typedef typename RemoveRef<MC>::Type MatrixC;
    typedef typename MatrixC::IndexType  IndexType;

#   ifndef NDEBUG
    int kA = (noTransA) ? A.numCols() : A.numRows();
    int kB = (noTransB) ? B.numRows() : B.numCols();
    ASSERT(kA==kB);
#   endif

    IndexType m = (noTransA) ? A.numRows() : A.numCols();
    IndexType n = (noTransB) ? B.numCols() : B.numRows();

    if (C.order()!=A.order()) {
        transposeA = Transpose(transposeA ^ Trans);
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

    if(DEBUGCLOSURE::identical(A, C)) {

        ASSERT( beta==BETA(0) );
        ASSERT( transposeB==NoTrans );

        FLENS_BLASLOG_SETTAG("--> ");
        FLENS_BLASLOG_BEGIN_TBMM(side, transA, alpha, B, C);

#       ifdef HAVE_CXXBLAS_TBMM
        cxxblas::tbmm(C.order(), Right,
                      Upper, transposeA, NonUnit,
                      C.numRows(), C.numCols(),
                      IndexType(0),
                      alpha,
                      B.data(), B.leadingDimension(),
                      C.data(), C.leadingDimension());
#       else
        ASSERT(0);
#       endif

        FLENS_BLASLOG_END;
        FLENS_BLASLOG_UNSETTAG;

        return;
    }

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_GBMM(transposeA, transposeB, alpha, A, B, beta, C);

#   ifdef HAVE_CXXBLAS_GBMM
    cxxblas::gbmm(C.order(), cxxblas::Side::Right,
                  transposeB, transposeA,
                  A.numRows(), A.numCols(),
                  IndexType(0), IndexType(0),
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


//-- gbmm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
typename RestrictTo<IsGbMatrix<MA>::value
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
        C.resize(m, n);
    }

    ASSERT(!DEBUGCLOSURE::identical(B, C));

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

//-- gbmm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsGbMatrix<MB>::value
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
        C.resize(m, n);
    }

    ASSERT(!DEBUGCLOSURE::identical(A, C));

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
        typename MatrixC::NoView _C;
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
    cxxblas::gemm(C.order(),
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

//== HermitianMatrix - GeneralMatrix products ==================================

//-- hbmm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
typename RestrictTo<IsHbMatrix<MA>::value
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
#   ifndef NDEBUG
    ASSERT(MC::order==MB::order);
    if (side==Left) {
        ASSERT(A.dim()==B.numRows());
    } else {
        ASSERT(B.numCols()==A.dim());
    }
#   endif

    StorageUpLo upLo = (C.order()==A.order())
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

    ASSERT(!DEBUGCLOSURE::identical(B, C));

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
    ASSERT(C.order()==B.order());
    if (side==Left) {
        ASSERT(A.dim()==B.numRows());
    } else {
        ASSERT(B.numCols()==A.dim());
    }
#   endif

    StorageUpLo upLo = (C.order()==A.order())
                     ? A.upLo()
                     : StorageUpLo(! A.upLo());

    IndexType m = (side==Left) ? A.dim() : B.numRows();
    IndexType n = (side==Left) ? B.numCols() : A.dim();

    ASSERT((beta==static_cast<BETA>(0)) || (C.numRows()==m));
    ASSERT((beta==static_cast<BETA>(0)) || (C.numCols()==n));

    if ((C.numRows()!=m) || (C.numCols()!=n)) {
        ASSERT(C.numRows()==0 && C.numCols()==0);
        C.resize(m,n);
    }

#   ifdef HAVE_CXXBLAS_HEMM
    cxxblas::hemm(C.order(), side,
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

//== SymmetricMatrix - GeneralMatrix products ==================================

//-- sbmm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
typename RestrictTo<IsSbMatrix<MA>::value
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
#   ifndef NDEBUG
    ASSERT(MC::order==MB::order);
    if (side==Left) {
        ASSERT(A.dim()==B.numRows());
    } else {
        ASSERT(B.numCols()==A.dim());
    }
#   endif

    StorageUpLo upLo = (C.order()==A.order())
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

    ASSERT(!DEBUGCLOSURE::identical(B, C));

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

//-- symm
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
    ASSERT(C.order()==B.order());
    if (side==Left) {
        ASSERT(A.dim()==B.numRows());
    } else {
        ASSERT(B.numCols()==A.dim());
    }
#   endif

    StorageUpLo upLo = (C.order()==A.order())
                     ? A.upLo()
                     : StorageUpLo(! A.upLo());

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
        typename MatrixC::NoView _C;
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
    cxxblas::symm(C.order(), side,
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


//== TriangularMatrix - GeneralMatrix products =================================

//-- tbmm
template <typename ALPHA, typename MA, typename MB>
    typename RestrictTo<IsTbMatrix<MA>::value
                     && IsGeMatrix<MB>::value,
             void>::Type
    mm(Side             side,
       Transpose        transA,
       const ALPHA      &alpha,
       const MA         &A,
       MB               &&B)
{
#   ifndef NDEBUG
    ASSERT(B.order()==A.order());
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

//-- trmm
template <typename ALPHA, typename MA, typename MB>
    typename RestrictTo<IsTrMatrix<MA>::value
                     && IsGeMatrix<MB>::value,
             void>::Type
    mm(Side             side,
       Transpose        transA,
       const ALPHA      &alpha,
       const MA         &A,
       MB               &&B)
{
#   ifndef NDEBUG
    ASSERT(B.order()==A.order());
    if (side==Left) {
        assert(A.dim()==B.numRows());
    } else {
        assert(B.numCols()==A.dim());
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_TRMM(side, transA, alpha, A, B);

#   ifdef HAVE_CXXBLAS_TRMM
    cxxblas::trmm(B.order(), side,
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

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL3_MM_TCC
