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

#ifndef FLENS_BLAS_LEVEL1_COPY_TCC
#define FLENS_BLAS_LEVEL1_COPY_TCC 1

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

//-- copy
template <typename VX, typename VY>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
copy(const VX &x, VY &&y)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_COPY(x, y);

//
//  Resize left hand size if needed.  This is *usually* only alloweded
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

#   ifdef HAVE_CXXBLAS_COPY
    cxxblas::copy(x.length(), x.data(), x.stride(), y.data(), y.stride());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- BLAS Level 1 extensions ---------------------------------------------------

//-- gecopy
template <typename MA, typename MB>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const MA &A, MB &&B)
{
    typedef typename RemoveRef<MA>::Type MatrixA;
    typedef typename RemoveRef<MB>::Type MatrixB;

//
//  check if this is an inplace transpose of A
//
    if (trans==Trans || trans==ConjTrans) {
        if (DEBUGCLOSURE::identical(A, B)) {
#           ifndef FLENS_DEBUG_CLOSURES
//
//          temporaries are not allowed
//
            ASSERT(A.numRows()==A.numCols());
            cxxblas::gecotr(B.order(), trans, B.numRows(), B.numCols(),
                            B.data(), B.leadingDimension());
            return;
#           else
//
//          temporaries are allowed: check if this requires a temporary
//
            if (A.numRows()!=A.numCols()) {
                typedef typename RemoveRef<MA>::Type   MatrixA;

                typename Result<MatrixA>::Type _A = A;
                FLENS_BLASLOG_TMP_ADD(_A);

                copy(trans, _A, B);

                FLENS_BLASLOG_TMP_REMOVE(_A, A);
                return;
            } else {
//
//              otherwise perform inplace transpose
//
                cxxblas::gecotr(B.order(), trans, B.numRows(), B.numCols(),
                                B.data(), B.leadingDimension());
                return;
            }
#           endif
        }
    }

//
//  Resize left hand size if needed.  This is *usually* only alloweded
//  when the left hand side is an empty matrix (such that it is no actual
//  resizing but rather an initialization).
//
    if (trans==NoTrans) {
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
            ASSERT(B.numRows()==0 && A.numCols()==0);
#           else
            if (B.numRows()!=0 || B.numCols()!=0) {
                FLENS_BLASLOG_RESIZE_MATRIX(B, A.numCols(), A.numRows());
            }
#           endif
            B.resize(A.numCols(), A.numRows(),
                     A.firstCol(), A.firstRow());
        }
    }

    trans = (A.order()==B.order())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MCOPY(trans, A, B);

#   ifdef HAVE_CXXBLAS_GECOPY
    cxxblas::gecopy(B.order(), trans, B.numRows(), B.numCols(),
                    A.data(), A.leadingDimension(),
                    B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- gbcopy
template <typename MA, typename MB>
typename RestrictTo<IsGbMatrix<MA>::value
                 && IsGbMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const MA &A, MB &&B)
{
    typename GbMatrix<MB>::ElementType  Zero(0);
//
//  check if this is an inplace transpose of A
//
    if (trans==Trans || trans==ConjTrans) {
        if (DEBUGCLOSURE::identical(A, B)) {
#           ifndef FLENS_DEBUG_CLOSURES
//
//          temporaries are not allowed
//
            ASSERT(A.numRows()==A.numCols());
            ASSERT(A.numSubDiags()==A.numSuperDiags());

            gbcotr(MB::order, trans, B.numRows(), B.numCols(),
                   B.numSubDiags(), B.numSuperDiags(),
                   B.data(), B.leadingDimension());
            return;
#           else
//
//          temporaries are allowed: check if this requires a temporary
//
            if ((A.numRows()!=A.numCols()) || 
                (A.numSubDiags()!=A.numSuperDiags())) {
                typename Result<GbMatrix<MA> >::Type _A = A;
                FLENS_BLASLOG_TMP_ADD(_A);

                copy(trans, _A, B);

                FLENS_BLASLOG_TMP_REMOVE(_A, A);
                return;
            } else {
//
//              otherwise perform inplace transpose
//
                gbcotr(MB::order, trans, B.numRows(), B.numCols(),
                       B.numSubDiags(), B.numSuperDiags(),
                       B.data(), B.leadingDimension());
                return;
            }
#           endif
        }
    }

//
//  Resize left hand size if needed.  This is *usually* only alloweded
//  when the left hand side is an empty matrix (such that it is no actual
//  resizing but rather an initialization).
//
    if ((trans==NoTrans) || (trans==Conj)) {
        if ((A.numRows()!=B.numRows()) || (A.numCols()!=B.numCols()) || 
            (A.numSubDiags()>B.numSubDiags()) || 
            (A.numSuperDiags()>B.numSuperDiags())) {
#           ifndef FLENS_DEBUG_CLOSURES
            ASSERT(B.numRows()==0 || B.numCols()==0);
#           else
            if (B.numRows()!=0 || B.numCols()!=0) {
                FLENS_BLASLOG_RESIZE_GBMATRIX(B, 
                                              A.numRows(), A.numCols(),
                                              A.numSubDiags(), 
                                              A.numSuperDiags());
            }
#           endif
            B.resize(A);
        }
    } else {
        if ((A.numRows()!=B.numCols())  || (A.numCols()!=B.numRows()) ||
            (A.numSubDiags()>B.numSuperDiags()) || 
            (A.numSuperDiags()>B.numSubDiags())) {
#           ifndef FLENS_DEBUG_CLOSURES
            ASSERT(B.numRows()==0 || B.numCols()==0);
#           else
            if (B.numRows()!=0 || B.numCols()!=0) {
                FLENS_BLASLOG_RESIZE_GBMATRIX(B, 
                                              A.numCols(), A.numRows(), 
                                              A.numSuperDiags(), 
                                              A.numSubDiags());
            }
#           endif

            B.resize(A.numCols(), A.numRows(),
                     A.numSuperDiags(), A.numSubDiags(),
                     A.firstIndex());
        }
    }

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MCOPY(trans, A, B);

    typedef typename GbMatrix<MB>::IndexType  IndexType;
    
    const IndexType numSubDiags   = ((trans==NoTrans) || (trans==Conj)) 
                                     ? A.numSubDiags() : A.numSuperDiags();
    const IndexType numSuperDiags = ((trans==NoTrans) || (trans==Conj)) 
                                     ? A.numSuperDiags() : A.numSubDiags();
    const IndexType Bshift = (MB::order==RowMajor) 
                               ? B.numSubDiags() - numSubDiags
                               : B.numSuperDiags() - numSuperDiags;
                             
    trans = (A.order()==B.order())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);
          
#   ifdef HAVE_CXXBLAS_GECOPY
    cxxblas::gbcopy(MB::order, trans,
                    B.numRows(), B.numCols(),
                    numSubDiags, numSuperDiags,
                    A.data(), A.leadingDimension(),
                    B.data()+Bshift, B.leadingDimension());
#   else
    ASSERT(0);
#   endif
   
    for(IndexType i = -B.numSubDiags(); i < -numSubDiags; ++i)
        B.viewDiag(i) = Zero;
    for(IndexType i = numSuperDiags+1; i <= B.numSuperDiags(); ++i)
        B.viewDiag(i) = Zero;

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- tbcopy
template <typename MA, typename MB>
typename RestrictTo<IsTbMatrix<MA>::value
                 && IsTbMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const MA &A, MB &&B)
{
    typename TbMatrix<MB>::ElementType  Zero(0), One(1);
    
    // Copy non-Unit diagonal into unit diagonal not possible
    ASSERT((A.diag()!=Unit) || (A.diag()!=B.diag()));  
  
//
//  Resize left hand size if needed.  This is *usually* only alloweded
//  when the left hand side is an empty matrix (such that it is no actual
//  resizing but rather an initialization).
//
    if ((trans==NoTrans) || (trans==Conj)) {
        ASSERT((A.upLo()==B.upLo()) || (B.dim()==0));
        if ((A.dim()!=B.dim())  || 
            (A.numOffDiags()>B.numOffDiags()) ) {
#           ifndef FLENS_DEBUG_CLOSURES
            ASSERT(B.dim()==0);
#           else
            if (B.dim()!=0 ) {
                FLENS_BLASLOG_RESIZE_TBMATRIX(B, 
                                              A.dim(), 
                                              A.numOffDiags());
            }
#           endif
            B.resize(A);

        }
    } else {
        ASSERT((A.upLo()!=B.upLo()) || (B.dim()==0));
        if ((A.dim()!=B.dim()) ||
            (A.numOffDiags()>B.numOffDiags())) {
#           ifndef FLENS_DEBUG_CLOSURES
            ASSERT(B.dim()==0);
#           else
            if (B.dim()!=0) {
                FLENS_BLASLOG_RESIZE_TBMATRIX(B, 
                                              A.dim(), 
                                              A.numOffDiags());
            }
#           endif

            B.resize(A.dim(), 
                     A.numOffDiags(), 
                     A.firstIndex());
        }
    }


    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MCOPY(trans, A, B);


    
    typedef typename TbMatrix<MB>::IndexType  IndexType;
    
    const IndexType numSubDiags   = ((trans==NoTrans) || (trans==Conj)) 
                                     ? A.numSubDiags() : A.numSuperDiags();
    const IndexType numSuperDiags = ((trans==NoTrans) || (trans==Conj)) 
                                     ? A.numSuperDiags() : A.numSubDiags();

    const IndexType Bshift = (MB::order==RowMajor) 
                               ? B.numSubDiags() - numSubDiags
                               : B.numSuperDiags() - numSuperDiags;

                               
    trans = (A.order()==B.order())
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);
          
#   ifdef HAVE_CXXBLAS_GECOPY
          
    cxxblas::gbcopy(MB::order, trans,
                    B.numRows(), B.numCols(),
                    numSubDiags, numSuperDiags,
                    A.data(), A.leadingDimension(),
                    B.data()+Bshift, B.leadingDimension());
#   else
    ASSERT(0);
#   endif
   
    if ((A.diag()==Unit) && (B.diag()==NonUnit))
      B.viewDiag(0) = One;
    
    for(IndexType i = -B.numSubDiags(); i < -numSubDiags; ++i)
        B.viewDiag(i) = Zero;
    for(IndexType i = numSuperDiags+1; i <= B.numSuperDiags(); ++i)
        B.viewDiag(i) = Zero;
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- tpcopy
template <typename MA, typename MB>
typename RestrictTo<IsTpMatrix<MA>::value
                 && IsTpMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const MA &A, MB &&B)
{
    ASSERT(A.diag()==B.diag()); 
    ASSERT(((A.upLo()==B.upLo()) && ((trans==NoTrans) || (trans==Conj))) || 
           ((A.upLo()!=B.upLo()) && ((trans==Trans) || (trans==ConjTrans))) );
//
//  Resize left hand size if needed.  This is *usually* only alloweded
//  when the left hand side is an empty matrix (such that it is no actual
//  resizing but rather an initialization).
//
    if (A.dim()!=B.dim()) {
#           ifndef FLENS_DEBUG_CLOSURES
            ASSERT(B.dim()==0);
#           else
            if (B.dim()!=0) {
                FLENS_BLASLOG_RESIZE_MATRIX(B, A.dim(), A.dim());
            }
#           endif
            B.resize(A);
    }
    
   trans = (A.order()==B.order())
              ? Transpose(trans ^ NoTrans)
              : Transpose(trans ^ Trans);
              
#   ifdef HAVE_CXXBLAS_TPCOPY
    cxxblas::tpcopy(MB::order, B.upLo(), trans, B.diag(),
                   B.dim(), A.data(), B.data());
 #   else
    ASSERT(0);
#   endif
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;

   return;
}

//-- sbcopy
template <typename MA, typename MB>
typename RestrictTo<IsSbMatrix<MA>::value
                 && IsSbMatrix<MB>::value,
         void>::Type
copy(const MA &A, MB &&B)
{
    typename SbMatrix<MB>::ElementType  Zero(0);
    
    Transpose trans = (A.upLo()==B.upLo()) ? NoTrans : Trans;

//
//  Resize left hand size if needed.  This is *usually* only alloweded
//  when the left hand side is an empty matrix (such that it is no actual
//  resizing but rather an initialization).
//
    if ((B.dim()!=A.dim()) || (B.numOffDiags()<A.numOffDiags())) {
#       ifndef FLENS_DEBUG_CLOSURES
        ASSERT(B.dim()==0);
#       else
        if (B.dim()!=0) {
            FLENS_BLASLOG_RESIZE_SBMATRIX(B, A.dim(), A.numOffDiags());
        }
#       endif
        B.resize(A);
    }

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_COPY(A, B);

    typedef typename SbMatrix<MB>::IndexType  IndexType;

    const IndexType numSubDiags   = (trans==NoTrans) 
                                     ? ((A.upLo()==Lower) ? A.numOffDiags() : 0)
                                     : ((A.upLo()==Upper) ? A.numOffDiags() : 0);
    const IndexType numSuperDiags = (trans==NoTrans)
                                     ? ((A.upLo()==Upper) ? A.numOffDiags() : 0) 
                                     : ((A.upLo()==Lower) ? A.numOffDiags() : 0);
    
    const IndexType Bshift = (MB::order==RowMajor) 
                               ? ((B.upLo()==Lower) ? B.numOffDiags() : 0) - numSubDiags
                               : ((B.upLo()==Upper) ? B.numOffDiags() : 0) - numSuperDiags; 
                               
    trans = (A.order()==B.order())
              ? Transpose(trans ^ NoTrans)
              : Transpose(trans ^ Trans);
              
#   ifdef HAVE_CXXBLAS_GBCOPY              
    cxxblas::gbcopy(MB::order, trans,
                    B.dim(), B.dim(),
                    numSubDiags, numSuperDiags,
                    A.data(), A.leadingDimension(),
                    B.data()+Bshift, B.leadingDimension());   
#   else
    ASSERT(0);
#   endif   
    
    for(IndexType i = A.numOffDiags()+1; i <= B.numOffDiags(); ++i)
            B.viewDiag(i) = Zero;

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- hbcopy
template <typename MA, typename MB>
typename RestrictTo<IsHbMatrix<MA>::value
                 && IsHbMatrix<MB>::value,
         void>::Type
copy(const MA &A, MB &&B)
{
    typename HbMatrix<MB>::ElementType  Zero(0);

    Transpose trans = (A.upLo()==B.upLo()) ? NoTrans : ConjTrans;
    
//
//  Resize left hand size if needed.  This is *usually* only alloweded
//  when the left hand side is an empty matrix (such that it is no actual
//  resizing but rather an initialization).
//
    if ((B.dim()!=A.dim()) || (B.numOffDiags()<A.numOffDiags())) {
#       ifndef FLENS_DEBUG_CLOSURES
        ASSERT(B.dim()==0);
#       else
        if (B.dim()!=0) {
            FLENS_BLASLOG_RESIZE_HBMATRIX(B, A.dim(), A.numOffDiags());
        }
#       endif
        B.resize(A);
    }

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_COPY(A, B);

 
    typedef typename SbMatrix<MB>::IndexType  IndexType;

    const IndexType numSubDiags   = (trans==NoTrans) 
                                     ? ((A.upLo()==Lower) ? A.numOffDiags() : 0)
                                     : ((A.upLo()==Upper) ? A.numOffDiags() : 0);
    const IndexType numSuperDiags = (trans==NoTrans)
                                     ? ((A.upLo()==Upper) ? A.numOffDiags() : 0) 
                                     : ((A.upLo()==Lower) ? A.numOffDiags() : 0);
    
    const IndexType Bshift = (MB::order==RowMajor) 
                               ? ((B.upLo()==Lower) ? B.numOffDiags() : 0) - numSubDiags
                               : ((B.upLo()==Upper) ? B.numOffDiags() : 0) - numSuperDiags;  
                               
    trans = (A.order()==B.order())
              ? Transpose(trans ^ NoTrans)
              : Transpose(trans ^ Trans);
#   ifdef HAVE_CXXBLAS_GBCOPY                               
    cxxblas::gbcopy(MB::order, trans,
                    B.dim(), B.dim(),
                    numSubDiags, numSuperDiags,
                    A.data(), A.leadingDimension(),
                    B.data()+Bshift, B.leadingDimension());       
#   else
    ASSERT(0);
#   endif
       
    for(IndexType i = A.numOffDiags()+1; i <= B.numOffDiags(); ++i)
            B.viewDiag(i) = Zero;
    
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- hpcopy
template <typename MA, typename MB>
typename RestrictTo<IsHpMatrix<MA>::value
                 && IsHpMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const MA &A, MB &&B)
{
    ASSERT(((A.upLo()==B.upLo()) && ((trans==NoTrans) || (trans==Conj))) || 
           ((A.upLo()!=B.upLo()) && ((trans==Trans) || (trans==ConjTrans))) );
//
//  Resize left hand size if needed.  This is *usually* only alloweded
//  when the left hand side is an empty matrix (such that it is no actual
//  resizing but rather an initialization).
//
    if (A.dim()!=B.dim()) {
#           ifndef FLENS_DEBUG_CLOSURES
            ASSERT(B.dim()==0);
#           else
            if (B.dim()!=0) {
                FLENS_BLASLOG_RESIZE_MATRIX(B, A.dim(), A.dim());
            }
#           endif
            B.resize(A);
    }
    
   trans = (A.order()==B.order())
              ? Transpose(trans ^ NoTrans)
              : Transpose(trans ^ Trans);

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_COPY(A, B);
   
#   ifdef HAVE_CXXBLAS_TPCOPY
    cxxblas::tpcopy(MB::order, B.upLo(), trans, B.diag(),
                   B.dim(), A.data(), B.data());
#   else
    ASSERT(0);
#   endif
    
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
    return;
}

//-- spcopy
template <typename MA, typename MB>
typename RestrictTo<IsSpMatrix<MA>::value
                 && IsSpMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const MA &A, MB &&B)
{
    ASSERT(((A.upLo()==B.upLo()) && ((trans==NoTrans) || (trans==Conj))) || 
           ((A.upLo()!=B.upLo()) && ((trans==Trans) || (trans==ConjTrans))) );
//
//  Resize left hand size if needed.  This is *usually* only alloweded
//  when the left hand side is an empty matrix (such that it is no actual
//  resizing but rather an initialization).
//
    if (A.dim()!=B.dim()) {
#           ifndef FLENS_DEBUG_CLOSURES
            ASSERT(B.dim()==0);
#           else
            if (B.dim()!=0) {
                FLENS_BLASLOG_RESIZE_MATRIX(B, A.dim(), A.dim());
            }
#           endif
            B.resize(A);
    }
    
   trans = (A.order()==B.order())
              ? Transpose(trans ^ NoTrans)
              : Transpose(trans ^ Trans);
    
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_COPY(A, B);
    
#   ifdef HAVE_CXXBLAS_TPCOPY
    cxxblas::tpcopy(MB::order, B.upLo(), trans, B.diag(),
                   B.dim(), A.data(), B.data());
#   else
    ASSERT(0);
#   endif
    
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
    return;
}


//-- trcopy
template <typename MA, typename MB>
typename RestrictTo<IsTrMatrix<MA>::value
                 && IsTrMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const MA &A, MB &&B)
{
//
//  Resize left hand size if needed.  This is *usually* only alloweded
//  when the left hand side is an empty matrix (such that it is no actual
//  resizing but rather an initialization).
//
    if (trans==NoTrans) {
        if ((A.numRows()!=B.numRows()) || (A.numCols()!=B.numCols())) {
#           ifndef FLENS_DEBUG_CLOSURES
            ASSERT(B.numRows()==0 || B.numCols()==0);
#           else
            if (B.numRows()!=0 && B.numCols()!=0) {
                FLENS_BLASLOG_RESIZE_MATRIX(B, A.numRows(), A.numCols());
            }
#           endif
            B.resize(A);
        }
    } else {
        if ((A.numRows()!=B.numCols()) || (A.numCols()!=B.numRows())) {
#           ifndef FLENS_DEBUG_CLOSURES
            ASSERT(B.numRows()==0 || B.numCols()==0);
#           else
            if (B.numRows()!=0 && B.numCols()!=0) {
                FLENS_BLASLOG_RESIZE_MATRIX(B, A.numCols(), A.numRows());
            }
#           endif
            B.resize(A.numCols(), A.numRows(),
                     A.firstCol(), A.firstRow());
        }
    }

#   ifndef NDEBUG
    if (trans==NoTrans || trans==Conj) {
        ASSERT(A.upLo()==B.upLo());
    } else {
        ASSERT(A.upLo()!=B.upLo());
    }
#   endif

    // TODO: make this assertion unnecessary
    ASSERT(A.order()==B.order());
    ASSERT(A.diag()==B.diag());

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_MCOPY(trans, A, B);

#   ifdef HAVE_CXXBLAS_TRCOPY
    cxxblas::trcopy(B.order(), B.upLo(), trans, B.diag(),
                    B.numRows(), B.numCols(),
                    A.data(), A.leadingDimension(),
                    B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- sycopy
template <typename MA, typename MB>
typename RestrictTo<IsSyMatrix<MA>::value
                 && IsSyMatrix<MB>::value,
         void>::Type
copy(const MA &A, MB &&B)
{
//
//  Resize left hand size if needed.  This is *usually* only alloweded
//  when the left hand side is an empty matrix (such that it is no actual
//  resizing but rather an initialization).
//
    if (B.dim()!=A.dim()) {
#       ifndef FLENS_DEBUG_CLOSURES
        ASSERT(B.dim()==0);
#       else
        if (B.dim()!=0) {
            FLENS_BLASLOG_RESIZE_MATRIX(B, A.dim(), A.dim());
        }
#       endif
        B.resize(A);
        B.upLo() = A.upLo();
    }

    ASSERT(A.upLo()==B.upLo());
    // TODO: make this assertion unnecessary
    ASSERT(A.order()==B.order());

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_COPY(A, B);

#   ifdef HAVE_CXXBLAS_TRCOPY
    cxxblas::trcopy(B.order(), B.upLo(), NoTrans, NonUnit, B.dim(), B.dim(),
                    A.data(), A.leadingDimension(),
                    B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif
    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- Sparse BLAS extensions ----------------------------------------------------

//-- copy: GeCoordMatrix -> GeCCSMatrix
template <typename MA, typename MB>
typename RestrictTo<IsGeCoordMatrix<MA>::value
                 && IsGeCCSMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const MA &A, MB &&B)
{
    B.engine() = A.engine();
}

//-- copy: GeCoordMatrix -> GeCRSMatrix
template <typename MA, typename MB>
typename RestrictTo<IsGeCoordMatrix<MA>::value
                 && IsGeCRSMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const MA &A, MB &&B)
{
    B.engine() = A.engine();
}

//-- copy: SyCoordMatrix -> SyCCSMatrix
template <typename MA, typename MB>
typename RestrictTo<IsSyCoordMatrix<MA>::value
                 && IsSyCCSMatrix<MB>::value,
         void>::Type
copy(const MA &A, MB &&B)
{
    B.engine() = A.engine();
    B.upLo() = A.upLo();
}

//-- copy: SyCoordMatrix -> SyCRSMatrix
template <typename MA, typename MB>
typename RestrictTo<IsSyCoordMatrix<MA>::value
                 && IsSyCRSMatrix<MB>::value,
         void>::Type
copy(const MA &A, MB &&B)
{
    B.engine() = A.engine();
    B.upLo() = A.upLo();
}


//-- convenience extensions ----------------------------------------------------


//-- copy: HbMatrix -> GbMatrix
template <typename MA, typename MB>
typename RestrictTo<IsHbMatrix<MA>::value
                 && IsGbMatrix<MB>::value,
         void>::Type
copy(const MA &A, MB &&B)
{
    
    if ((A.numRows()!=B.numRows()) ||
        (A.numCols()!=B.numCols()) || 
        (A.numOffDiags() > B.numSubDiags()) || 
        (A.numOffDiags() > B.numSuperDiags())) {
#       ifndef FLENS_DEBUG_CLOSURES
        ASSERT(B.numRows()==0 || B.numCols()==0);
#       else
        if (B.numRows()!=0 && B.numCols()!=0) {
            FLENS_BLASLOG_RESIZE_GBMATRIX(B, A.numRows(), A.numCols(),
                                      A.numOffDiags(), A.numOffDiags());
        }
#       endif
        B.resize(A.numRows(), A.numCols(),
                 A.numOffDiags(), A.numOffDiags(),
                 A.firstIndex(), A.firstIndex());
    }

    if (A.upLo()==Upper) {
        B.upper() = A.triangular();
        B.strictLower() = conjgateTranspose(A.general().strictUpper());
    } else {
        B.lower() = A.triangular();
        B.strictUpper() = conjugateTranspose(A.general().strictLower());
    }
}


//-- copy: TbMatrix -> GbMatrix
template <typename MA, typename MB>
typename RestrictTo<IsTbMatrix<MA>::value
                 && IsGbMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const MA &A, MB &&B)
{
  
    typename GeMatrix<MB>::ElementType  Zero(0);
    typename GeMatrix<MB>::ElementType  One(1);

    if (trans==NoTrans) {
        if ((A.numRows()!=B.numRows()) 
            || (A.numCols()!=B.numCols())
            || (A.numSubDiags()>B.numSubDiags())
            || (A.numSuperDiags()>B.numSuperDiags())) {
#           ifndef FLENS_DEBUG_CLOSURES
            ASSERT((B.numRows()==0) && (B.numCols()==0));
#           else
            if (B.numRows()!=0 && B.numCols()!=0) {
                FLENS_BLASLOG_RESIZE_GBMATRIX(B, A.numRows(), A.numCols(),
                                A.numSubDiags(), A.numSuperDiags());
            }
#           endif
            B.resize(A.numRows(), A.numCols(),
                     A.numSubDiags(), A.numSuperDiags(),
                     A.firstIndex());
        }
    } else {
        if ((A.numRows()!=B.numCols()) 
            || (A.numCols()!=B.numRows())
            || (A.numSubDiags()>B.numSuperDiags())
            || (A.numSuperDiags()>B.numSubDiags())) {
#           ifndef FLENS_DEBUG_CLOSURES
            ASSERT((B.numRows()==0) && (B.numCols()==0));
#           else
            if (B.numRows()!=0 && B.numCols()!=0) {
                FLENS_BLASLOG_RESIZE_GBMATRIX(B, A.numCols(), A.numRows(),
                                A.numSuperDiags(), A.numSubDiags());
            }
#           endif
            B.resize(A.numCols(), A.numRows(),
                     A.numSuperDiags(), A.numSubDiags());
         }
    }
    typedef typename TbMatrix<MB>::IndexType  IndexType;
    if (trans==NoTrans) {
        if (A.upLo()==Upper) {
            B.upper() = A;
            B.strictLower() = Zero;
        } else {
            B.lower() = A;
            B.strictUpper() = Zero;
        }
        if (A.diag()==Unit)
            B.viewDiag(0)=One;
        
    } else if (trans==Conj) {
        if (A.upLo()==Upper) {
            B.upper() = conjugate(A);
            B.strictLower() = Zero;
        } else {
            B.lower() = conjugate(A);
            B.strictUpper() = Zero;
        }
        if (A.diag()==Unit)
            B.viewDiag(0)=One;
        
    } else if (trans==Trans) {
        if (A.upLo()==Upper) {
            B.lower() = transpose(A);
            B.strictUpper() = Zero;
        } else {
            B.upper() = transpose(A);
            B.strictLower() = Zero;
        }
        if (A.diag()==Unit)
            B.viewDiag(0)=One;      
    } else {
        if (A.upLo()==Upper) {
            B.lower() = conjugateTranspose(A);
            B.strictUpper() = Zero;
        } else {
            B.upper() = conjugateTranspose(A);
            B.strictLower() = Zero;
        }
        if (A.diag()==Unit)
            B.viewDiag(0)=One;      
    }
}

//-- copy: TrMatrix -> GeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsTrMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const MA &A, MB &&B)
{
    typedef typename RemoveRef<MA>::Type MatrixB;

    typename MatrixB::ElementType  Zero(0), One(1);

    if (trans==NoTrans) {
        if (A.numRows()!=B.numRows() && A.numCols()!=B.numCols()) {
#           ifndef FLENS_DEBUG_CLOSURES
            ASSERT(B.numRows()==0 || B.numCols()==0);
#           else
            if (B.numRows()!=0 && B.numCols()!=0) {
                FLENS_BLASLOG_RESIZE_MATRIX(B, A.numRows(), A.numCols());
            }
#           endif
            B.resize(A.numRows(), A.numCols(),
                     A.firstRow(), A.firstCol());
        }
    } else {
        if (A.numRows()!=B.numCols() && A.numCols()!=B.numRows()) {
#           ifndef FLENS_DEBUG_CLOSURES
            ASSERT(B.numRows()==0 || B.numCols()==0);
#           else
            if (B.numRows()!=0 && B.numCols()!=0) {
                FLENS_BLASLOG_RESIZE_MATRIX(B, A.numCols(), A.numRows());
            }
#           endif
            B.resize(A.numCols(), A.numRows(),
                     A.firstCol(), A.firstRow());
         }
    }

    if (trans==NoTrans) {
        if (A.upLo()==Upper) {
            if (A.diag()!=Unit) {
                B.upper() = A;
            } else {
                B.upperUnit() = A;
                B.diag(0) = One;
            }
            B.strictLower() = Zero;
        } else {
            if (A.diag()!=Unit) {
                B.lower() = A;
            } else {
                B.lowerUnit() = A;
                B.diag(0) = One;
            }
            B.strictUpper() = Zero;
        }
    } else if (trans==Trans) {
        if (A.upLo()==Upper) {
            if (A.diag()!=Unit) {
                B.lower() = transpose(A);
            } else {
                B.lowerUnit() = transpose(A);
                B.diag(0) = One;
            }
            B.strictUpper() = Zero;
        } else {
            if (A.diag()!=Unit) {
                B.upper() = transpose(A);
            } else {
                B.upperUnit() = transpose(A);
                B.diag(0) = One;
            }
            B.strictLower() = Zero;
        }
    } else {
        ASSERT(0);
    }
}

//-- copy: SbMatrix -> GbMatrix
template <typename MA, typename MB>
typename RestrictTo<IsSbMatrix<MA>::value
                 && IsGbMatrix<MB>::value,
         void>::Type
copy(const MA &A, MB &&B)
{
    
    if ((A.numRows()!=B.numRows()) ||
        (A.numCols()!=B.numCols()) || 
        (A.numOffDiags() > B.numSubDiags()) || 
        (A.numOffDiags() > B.numSuperDiags())) {
#       ifndef FLENS_DEBUG_CLOSURES
        ASSERT(B.numRows()==0 || B.numCols()==0);
#       else
        if (B.numRows()!=0 && B.numCols()!=0) {
            FLENS_BLASLOG_RESIZE_GBMATRIX(B, A.numRows(), A.numCols(),
                                      A.numOffDiags(), A.numOffDiags());
        }
#       endif
        B.resize(A.numRows(), A.numCols(),
                 A.numOffDiags(), A.numOffDiags(),
                 A.firstIndex(), A.firstIndex());
    }

    if (A.upLo()==Upper) {
        B.upper() = A.triangular();
        B.strictLower() = transpose(A.triangular().general().strictUpper());
    } else {
        B.lower() = A.triangular();
        B.strictUpper() = transpose(A.triangular().general().strictLower());        

    }
}


//-- copy: SyMatrix -> GeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsSyMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
copy(const MA &A, MB &&B)
{
    if (A.numRows()!=B.numRows() && A.numCols()!=B.numCols()) {
#       ifndef FLENS_DEBUG_CLOSURES
        ASSERT(B.numRows()==0 || B.numCols()==0);
#       else
        if (B.numRows()!=0 && B.numCols()!=0) {
            FLENS_BLASLOG_RESIZE_MATRIX(B, A.numRows(), A.numCols());
        }
#       endif
        B.resize(A.numRows(), A.numCols(),
                 A.firstRow(), A.firstCol());
    }

    if (A.upLo()==Upper) {
        B.upper() = A.general().upper();
        B.strictLower() = transpose(A.general().strictUpper());
    } else {
        B.lower() = A.general().lower();
        B.strictUpper() = transpose(A.general().strictLower());
    }
}

//-- copy: GeCoordMatrix -> GeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsGeCoordMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const MA &A, MB &&B)
{
    typedef typename MA::ElementType  ElementType;

    B.resize(A.numRows(), A.numCols(),
             A.firstRow(), A.firstCol(),
             ElementType(0));

    const auto &coord = A.engine().coordVector();

    for (size_t k=0; k<coord.size(); ++k) {
        B(coord[k].row, coord[k].col) += coord[k].value;
    }
}

//-- copy: GeCCSMatrix -> GeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsGeCCSMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const MA &A, MB &&B)
{
    typedef typename MA::IndexType    IndexType;
    typedef typename MA::ElementType  ElementType;

    B.resize(A.numRows(), A.numCols(),
             A.firstRow(), A.firstCol(),
             ElementType(0));

    const auto &cols = A.engine().cols();
    const auto &rows = A.engine().rows();
    const auto &vals = A.engine().values();

    for (IndexType j=cols.firstIndex(); j<cols.lastIndex(); ++j) {
        for (IndexType k=cols(j); k<cols(j+1); ++k) {
            B(rows(k), j) = vals(k);
        }
    }
}

//-- copy: GeCRSMatrix -> GeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsGeCRSMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const MA &A, MB &&B)
{
    typedef typename MA::IndexType    IndexType;
    typedef typename MA::ElementType  ElementType;

    B.resize(A.numRows(), A.numCols(),
             A.firstRow(), A.firstCol(),
             ElementType(0));

    const auto &rows = A.engine().rows();
    const auto &cols = A.engine().cols();
    const auto &vals = A.engine().values();

    for (IndexType i=rows.firstIndex(); i<rows.lastIndex(); ++i) {
        for (IndexType k=rows(i); k<rows(i+1); ++k) {
            B(i,cols(k)) = vals(k);
        }
    }
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_COPY_TCC
