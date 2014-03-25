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

//-- BLAS Level 1 --------------------------------------------------------------

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
//      So we allow  y += alpha*x  for an empty vector y
//
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;
        const T  Zero(0);

        y.resize(x, Zero);
    }
    ASSERT(y.length()==x.length());

#   ifdef HAVE_CXXBLAS_AXPY
    cxxblas::axpy(x.length(), alpha,
                  x.data(), x.stride(),
                  y.data(), y.stride());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- axpy
template <typename ALPHA, typename VX, typename VY>
typename RestrictTo<IsTinyVector<VX>::value
                 && IsTinyVector<VY>::value,
         void>::Type
axpy(const ALPHA &alpha, const VX &x, VY &&y)
{
    typedef typename VX::ElementType       TX;
    typedef typename RemoveRef<VY>::Type   VectorY;
    typedef typename VectorY::ElementType  TY;

    const int n    = VX::Engine::length;
    const int incX = VX::Engine::stride;
    const int incY = VectorY::Engine::stride;

    cxxblas::axpy<n, ALPHA, TX, incX, TY, incY>(alpha, x.data(), y.data());
}

//-- BLAS Level 1 extensions ---------------------------------------------------

//== GeneralMatrix

//-- diagaxpy
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsDiagMatrix<MA>::value
                 && IsDiagMatrix<MB>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha, const MA &A, MB &&B)
{
    if ( trans==Conj || trans==ConjTrans) {
        B.diag() += alpha*conjugate(A.diag());
        return;
    }
    
    B.diag() += alpha*A.diag();
}

//-- gbaxpy
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsGbMatrix<MA>::value
                 && IsGbMatrix<MB>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha, const MA &A, MB &&B)
{
    typedef typename RemoveRef<MA>::Type   MatrixA;
    typedef typename MatrixA::IndexType  IndexType;
    
    if (B.numRows()==0 || B.numCols()==0) {
//
//      So we allow  B += alpha*A  for an empty matrix B
//
        if ((trans==NoTrans) || (trans==Conj)) {
            B.resize(A.numRows(), A.numCols(),
                     A.numSubDiags(), A.numSuperDiags(),
                     A.firstIndex());
        } else {
            B.resize(A.numCols(), A.numRows(),
                     A.numSuperDiags(), A.numSubDiags(),
                     A.firstIndex());
        }
    }

#   ifndef NDEBUG
    if ((trans==NoTrans) || (trans==Conj)) {
        ASSERT((A.numRows()==B.numRows()) && (A.numCols()==B.numCols()));
        ASSERT((A.numSubDiags()<=B.numSubDiags())
            && (A.numSuperDiags()<=B.numSuperDiags()));
    } else {
        ASSERT((A.numRows()==B.numCols()) && (A.numCols()==B.numRows()));
        ASSERT((A.numSubDiags()<=B.numSuperDiags())
            && (A.numSuperDiags()<=B.numSubDiags()));
    }
#   endif


#   ifndef NDEBUG
    if (trans==Trans || trans==ConjTrans) {
        ASSERT(!DEBUGCLOSURE::identical(A, B));
    }
#   endif


    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);

    const IndexType numSubDiags   = ((trans==NoTrans) || (trans==Conj))
                                  ? A.numSubDiags() : A.numSuperDiags();
    const IndexType numSuperDiags = ((trans==NoTrans) || (trans==Conj))
                                  ? A.numSuperDiags() : A.numSubDiags();
    const IndexType Bshift = (B.order()==RowMajor)
                           ? B.numSubDiags() - numSubDiags
                           : B.numSuperDiags() - numSuperDiags;


    trans = (A.order()==B.order())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);

#   ifdef HAVE_CXXBLAS_GBAXPY
    cxxblas::gbaxpy(B.order(), trans,
                    B.numRows(), B.numCols(),
                    numSubDiags, numSuperDiags, alpha,
                    A.data(), A.leadingDimension(),
                    B.data()+Bshift, B.leadingDimension());
#   else
    ASSERT(0);
#   endif
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- geaxpy
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha, const MA &A, MB &&B)
{
    typedef typename RemoveRef<MB>::Type   MatrixB;

    if (B.numRows()==0 || B.numCols()==0) {
//
//      So we allow  B += alpha*A  for an empty matrix B
//
        typedef typename MatrixB::ElementType  T;
        const T  Zero(0);

        if ((trans==NoTrans) || (trans==Conj)) {
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
//
//  If A and B are identical a temporary is needed if we want to use axpy
//  for B += alpha*A^T or B+= alpha*A^H
//
    typedef typename RemoveRef<MA>::Type   MatrixA;

    if ((trans==Trans || trans==ConjTrans) && DEBUGCLOSURE::identical(A, B)) {
        typename Result<MA>::Type _A = A;
        FLENS_BLASLOG_TMP_ADD(_A);

        copy(trans, A, _A);
        axpy(NoTrans, alpha, A, B);

        FLENS_BLASLOG_TMP_REMOVE(_A, A);
        return;
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);

#   ifdef HAVE_CXXBLAS_GEAXPY
    geaxpy(B.order(), trans, B.numRows(), B.numCols(), alpha,
           A.data(), A.leadingDimension(),
           B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- (tiny) geaxpy
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsGeTinyMatrix<MA>::value
                 && IsGeTinyMatrix<MB>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha, const MA &A, MB &&B)
{
    typedef typename RemoveRef<MA>::Type   MatrixA;
    typedef typename RemoveRef<MB>::Type   MatrixB;

    typedef typename MatrixA::ElementType  TA;
    typedef typename MatrixB::ElementType  TB;

    const int m   = MatrixA::Engine::numRows;
    const int n   = MatrixA::Engine::numCols;
    const int ldA = MatrixA::Engine::leadingDimension;
    const int ldB = MatrixB::Engine::leadingDimension;

    if (trans==Trans || trans==ConjTrans) {
        ASSERT(!DEBUGCLOSURE::identical(A, B));
    }

    cxxblas::geaxpy<m, n, ALPHA, TA, ldA, TB, ldB>(trans, alpha,
                                                   A.data(), B.data());
}

//== HermitianMatrix

//-- hbaxpy
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsHbMatrix<MA>::value
                 && IsHbMatrix<MB>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha, const MA &A, MB &&B)
{
    typedef typename RemoveRef<MA>::Type   MatrixA;
    typedef typename MatrixA::IndexType  IndexType;

    ASSERT(cxxblas::imag(alpha)==0);

    if (B.dim()==0) {
        B.resize(A);
    }

#   ifndef NDEBUG
    ASSERT(A.dim()==B.dim());
    ASSERT(A.numOffDiags()<=B.numOffDiags());
    if ((trans==NoTrans) || (trans==Conj)) {
        ASSERT(A.upLo() == B.upLo()) ;
    } else {
        ASSERT(A.upLo() != B.upLo()) ;
    }
#   endif


#   ifndef FLENS_DEBUG_CLOSURES
#   ifndef NDEBUG
    if (trans==Trans || trans==ConjTrans) {
        ASSERT(!DEBUGCLOSURE::identical(A, B));
    }
#   endif
#   else
//
//  If A and B are identical a temporary is needed if we want to use axpy
//  for B += alpha*A^T or B+= alpha*A^H
//
    if ((trans==Trans || trans==ConjTrans) && DEBUGCLOSURE::identical(A, B)) {
        typename SbMatrix<MA>::NoView _A;
        FLENS_BLASLOG_TMP_ADD(_A);

        copy(trans, A, _A);
        axpy(NoTrans, alpha, _A, B);

        FLENS_BLASLOG_TMP_REMOVE(_A, A);
        return;
    }
#   endif


    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);

    trans = (A.upLo() == B.upLo())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ ConjTrans);

    const IndexType numSubDiags   = (trans==NoTrans)
                                  ? ((A.upLo()==Lower) ? A.numOffDiags() : 0)
                                  : ((A.upLo()==Upper) ? A.numOffDiags() : 0);
    const IndexType numSuperDiags = (trans==NoTrans)
                                  ? ((A.upLo()==Upper) ? A.numOffDiags() : 0)
                                  : ((A.upLo()==Lower) ? A.numOffDiags() : 0);

    const IndexType Bshift = (B.order()==RowMajor)
                    ? ((B.upLo()==Lower) ? B.numOffDiags() : 0) - numSubDiags
                    : ((B.upLo()==Upper) ? B.numOffDiags() : 0) - numSuperDiags;

    trans = (A.order()==B.order())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);
#   ifdef HAVE_CXXBLAS_GBAXPY
    cxxblas::gbaxpy(B.order(), trans,
                    B.numRows(), B.numCols(),
                    numSubDiags, numSuperDiags, alpha,
                    A.data(), A.leadingDimension(),
                    B.data()+Bshift, B.leadingDimension());
#   else
    ASSERT(0);
#   endif
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//heaxpy
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsHeMatrix<MB>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha, const MA &A, MB &&B)
{
//
//  Resize left hand size if needed.  This is *usually* only alloweded
//  when the left hand side is an empty matrix (such that it is no actual
//  resizing but rather an initialization).
//
    if (B.dim()==0) {
        B.resize(A.dim());
    }

    ASSERT(A.dim()==B.dim());
    
    // alpha must be real
    ASSERT(cxxblas::imag(alpha)==0);
    
    trans = (A.upLo()==B.upLo())          
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ ConjTrans);
    
    trans = (A.order()==B.order())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);
          
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);

#   ifdef HAVE_CXXBLAS_TRAXPY
    cxxblas::traxpy(B.order(), B.upLo(), trans, NonUnit,
                     B.numRows(), B.numCols(), alpha,
                     A.data(), A.leadingDimension(),
                     B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- hpaxpy
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsHpMatrix<MA>::value
                 && IsHpMatrix<MB>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha, const MA &A, MB &&B)
{

    if (B.dim()==0) {
        B.resize(A);
    }

#   ifndef NDEBUG
    ASSERT(A.dim()==B.dim());
#   endif

#   ifndef NDEBUG
    if (trans==Trans || trans==ConjTrans) {
        ASSERT(!DEBUGCLOSURE::identical(A, B));
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);

    trans = (A.upLo() == B.upLo())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);

    trans = (A.order()==B.order())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);

#   ifdef HAVE_CXXBLAS_TPAXPY
    cxxblas::tpaxpy(B.order(), B.upLo(), trans,
                    NonUnit, B.dim(),
                    alpha, A.data(), B.data());
#   else
    ASSERT(0);
#   endif
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//== SymmetricMatrix

//-- sbaxpy
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsSbMatrix<MA>::value
                 && IsSbMatrix<MB>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha, const MA &A, MB &&B)
{
    typedef typename RemoveRef<MA>::Type   MatrixA;
    typedef typename MatrixA::IndexType  IndexType;

    if (B.dim()==0) {
        B.resize(A);
    }

#   ifndef NDEBUG
    ASSERT(A.dim()==B.dim());
    ASSERT(A.numOffDiags()<=B.numOffDiags());
    if ((trans==NoTrans) || (trans==Conj)) {
        ASSERT(A.upLo() == B.upLo()) ;
    } else {
        ASSERT(A.upLo() != B.upLo()) ;
    }
#   endif

#   ifndef FLENS_DEBUG_CLOSURES
#   ifndef NDEBUG
    if (trans==Trans || trans==ConjTrans) {
        ASSERT(!DEBUGCLOSURE::identical(A, B));
    }
#   endif
#   else
//
//  If A and B are identical a temporary is needed if we want to use axpy
//  for B += alpha*A^T or B+= alpha*A^H
//
    if ((trans==Trans || trans==ConjTrans) && DEBUGCLOSURE::identical(A, B)) {
        typename SbMatrix<MA>::NoView _A;
        FLENS_BLASLOG_TMP_ADD(_A);

        copy(trans, A, _A);
        axpy(NoTrans, alpha, _A, B);

        FLENS_BLASLOG_TMP_REMOVE(_A, A);
        return;
    }
#   endif


    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);

    trans = (A.upLo() == B.upLo())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);

    const IndexType numSubDiags = (trans==NoTrans)
                                ? ((A.upLo()==Lower) ? A.numOffDiags() : 0)
                                : ((A.upLo()==Upper) ? A.numOffDiags() : 0);
    const IndexType numSuperDiags = (trans==NoTrans)
                                  ? ((A.upLo()==Upper) ? A.numOffDiags() : 0)
                                  : ((A.upLo()==Lower) ? A.numOffDiags() : 0);

    const IndexType Bshift = (B.order()==RowMajor)
                    ? ((B.upLo()==Lower) ? B.numOffDiags() : 0) - numSubDiags
                    : ((B.upLo()==Upper) ? B.numOffDiags() : 0) - numSuperDiags;

    trans = (A.order()==B.order())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);

#   ifdef HAVE_CXXBLAS_GBAXPY
    cxxblas::gbaxpy(B.order(), trans,
                    B.numRows(), B.numCols(),
                    numSubDiags, numSuperDiags, alpha,
                    A.data(), A.leadingDimension(),
                    B.data()+Bshift, B.leadingDimension());
#   else
    ASSERT(0);
#   endif
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- spaxpy
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsSpMatrix<MA>::value
                 && IsSpMatrix<MB>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha, const MA &A, MB &&B)
{
    typedef typename RemoveRef<MA>::Type   MatrixA;
    typedef typename MatrixA::IndexType  IndexType;
    
    ASSERT(B.diag()==NonUnit);
    // TODO: Remove this condition
    ASSERT(A.diag()==NonUnit);

    if (B.dim()==0) {
        B.resize(A);
    }

#   ifndef NDEBUG
    ASSERT(A.dim()==B.dim());
    if ((trans==NoTrans) || (trans==Conj)) {
        ASSERT((A.numSubDiags()<=B.numSubDiags())
            && (A.numSuperDiags()<=B.numSuperDiags()));
    } else {
        ASSERT((A.numSubDiags()<=B.numSuperDiags())
            && (A.numSuperDiags()<=B.numSubDiags()));
    }
#   endif

#   ifndef NDEBUG
    if (trans==Trans || trans==ConjTrans) {
        ASSERT(!DEBUGCLOSURE::identical(A, B));
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);

    const IndexType numSubDiags = ((trans==NoTrans) || (trans==Conj))
                                ? A.numSubDiags() : A.numSuperDiags();
    const IndexType numSuperDiags = ((trans==NoTrans) || (trans==Conj))
                                  ? A.numSuperDiags() : A.numSubDiags();
    const IndexType Bshift = (B.order()==RowMajor)
                               ? B.numSubDiags() - numSubDiags
                               : B.numSuperDiags() - numSuperDiags;

    trans = (A.order()==B.order())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);

#   ifdef HAVE_CXXBLAS_GBAXPY
    cxxblas::gbaxpy(B.order(), trans,
                    B.numRows(), B.numCols(),
                    numSubDiags, numSuperDiags, alpha,
                    A.data(), A.leadingDimension(),
                    B.data()+Bshift, B.leadingDimension());
#   else
    ASSERT(0);
#   endif
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}


//-- syaxpy
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsSyMatrix<MA>::value
                 && IsSyMatrix<MB>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha, const MA &A, MB &&B)
{
//
//  Resize left hand size if needed.  This is *usually* only alloweded
//  when the left hand side is an empty matrix (such that it is no actual
//  resizing but rather an initialization).
//
    if (B.dim()==0) {
        B.resize(A.dim());
    }

    ASSERT(A.dim()==B.dim());
    trans = (A.upLo()==B.upLo())          
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);
    
    trans = (A.order()==B.order())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);
          
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);

#   ifdef HAVE_CXXBLAS_TRAXPY
    cxxblas::traxpy(B.order(), B.upLo(), trans, NonUnit,
                    B.numRows(), B.numCols(), alpha,
                    A.data(), A.leadingDimension(),
                    B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}


//== TriangularMatrix

//-- tbaxpy
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsTbMatrix<MA>::value
                 && IsTbMatrix<MB>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha, const MA &A, MB &&B)
{
    typedef typename RemoveRef<MA>::Type   MatrixA;
    typedef typename MatrixA::IndexType  IndexType;
    
    ASSERT(B.diag()==NonUnit);
    // TODO: Remove this condition
    ASSERT(A.diag()==NonUnit);

    if (B.dim()==0) {
        B.resize(A);
    }

#   ifndef NDEBUG
    ASSERT(A.dim()==B.dim());
    if ((trans==NoTrans) || (trans==Conj)) {
        ASSERT((A.numSubDiags()<=B.numSubDiags())
            && (A.numSuperDiags()<=B.numSuperDiags()));
    } else {
        ASSERT((A.numSubDiags()<=B.numSuperDiags())
            && (A.numSuperDiags()<=B.numSubDiags()));
    }
#   endif

#   ifndef NDEBUG
    if (trans==Trans || trans==ConjTrans) {
        ASSERT(!DEBUGCLOSURE::identical(A, B));
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);

    trans = (A.upLo() == B.upLo())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);

    const IndexType numSubDiags = (trans==NoTrans)
                                ? ((A.upLo()==Lower) ? A.numOffDiags() : 0)
                                : ((A.upLo()==Upper) ? A.numOffDiags() : 0);
    const IndexType numSuperDiags = (trans==NoTrans)
                                  ? ((A.upLo()==Upper) ? A.numOffDiags() : 0)
                                  : ((A.upLo()==Lower) ? A.numOffDiags() : 0);

    const IndexType Bshift = (B.order()==RowMajor)
                    ? ((B.upLo()==Lower) ? B.numOffDiags() : 0) - numSubDiags
                    : ((B.upLo()==Upper) ? B.numOffDiags() : 0) - numSuperDiags;

    trans = (A.order()==B.order())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);

#   ifdef HAVE_CXXBLAS_GBAXPY
    cxxblas::gbaxpy(B.order(), trans,
                    B.numRows(), B.numCols(),
                    numSubDiags, numSuperDiags, alpha,
                    A.data(), A.leadingDimension(),
                    B.data()+Bshift, B.leadingDimension());
#   else
    ASSERT(0);
#   endif
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- traxpy
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsTrMatrix<MA>::value
                 && IsTrMatrix<MB>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha, const MA &A, MB &&B)
{
    ASSERT(A.diag()==NonUnit);
    ASSERT(B.diag()==NonUnit);
    
#   ifndef NDEBUG
    if (A.upLo()==B.upLo()) {
        ASSERT(trans==NoTrans || trans==Conj) ;
    } else {
        ASSERT(trans==Trans || trans==ConjTrans) ;    
    }
#   endif

    trans = (A.order()==B.order())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);

#   ifdef HAVE_CXXBLAS_TRAXPY
    cxxblas::traxpy(B.order(), B.upLo(), trans, B.diag(),
                    B.numRows(), B.numCols(), alpha,
                    A.data(), A.leadingDimension(),
                    B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}



//-- tpaxpy
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsTpMatrix<MA>::value
                 && IsTpMatrix<MB>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha, const MA &A, MB &&B)
{
    ASSERT(B.diag()==NonUnit);
    // TODO: Remove this condition
    ASSERT(A.diag()==NonUnit);

    if (B.dim()==0) {
        B.resize(A);
    }

#   ifndef NDEBUG
    ASSERT(A.dim()==B.dim());
#   endif

#   ifndef NDEBUG
    if (trans==Trans || trans==ConjTrans) {
        ASSERT(!DEBUGCLOSURE::identical(A, B));
    }
#   endif

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MAXPY(trans, alpha, A, B);

    trans = (A.order()==B.order())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);

#   ifdef HAVE_CXXBLAS_TPAXPY
    cxxblas::tpaxpy(B.order(), B.upLo(), trans,
                    B.diag(), B.dim(),
                    alpha, A.data(), B.data());
#   else
    ASSERT(0);
#   endif
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_AXPY_TCC
