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
//  Resize left hand size if needed.  This is *usually* only allowed
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

//-- hecopy
template <typename MA, typename MB>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsHeMatrix<MB>::value,
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

//-- copy: HeCoordMatrix -> HeCCSMatrix
template <typename MA, typename MB>
typename RestrictTo<IsHeCoordMatrix<MA>::value
                 && IsHeCCSMatrix<MB>::value,
         void>::Type
copy(const MA &A, MB &&B)
{
    B.engine() = A.engine();
    B.upLo() = A.upLo();
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

//-- copy: GeCoordMatrix -> GeCRSMatrix
template <typename MA, typename MB>
typename RestrictTo<IsGeCoordMatrix<MA>::value
                 && IsGeCRSMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const MA &A, MB &&B)
{
    B.engine() = A.engine();
}

//-- copy: HeCoordMatrix -> HeCRSMatrix
template <typename MA, typename MB>
typename RestrictTo<IsHeCoordMatrix<MA>::value
                 && IsHeCRSMatrix<MB>::value,
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

//-- copy: HeCCSMatrix -> HeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsHeCCSMatrix<MA>::value
                 && IsHeMatrix<MB>::value,
         void>::Type
copy(const MA &A, MB &&B)
{
    typedef typename MA::IndexType    IndexType;
    typedef typename MA::ElementType  ElementType;

    if (B.dim()!=A.dim()) {
        ASSERT(B.dim()==0);
        B.resize(A.dim(), A.upLo(), A.indexBase(), ElementType(0));
        B.upLo() = A.upLo();
    }

    ASSERT(A.upLo()==B.upLo());

    const auto &cols = A.engine().cols();
    const auto &rows = A.engine().rows();
    const auto &vals = A.engine().values();

    for (IndexType j=cols.firstIndex(); j<cols.lastIndex(); ++j) {
        for (IndexType k=cols(j); k<cols(j+1); ++k) {
            B(rows(k), j) = vals(k);
        }
    }
}

//-- copy: SyCCSMatrix -> SyMatrix
template <typename MA, typename MB>
typename RestrictTo<IsSyCCSMatrix<MA>::value
                 && IsSyMatrix<MB>::value,
         void>::Type
copy(const MA &A, MB &&B)
{
    typedef typename MA::IndexType    IndexType;
    typedef typename MA::ElementType  ElementType;

    if (B.dim()!=A.dim()) {
        ASSERT(B.dim()==0);
        B.resize(A.dim(), A.upLo(), A.indexBase(), ElementType(0));
        B.upLo() = A.upLo();
    }

    ASSERT(A.upLo()==B.upLo());

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

//-- copy: HeCRSMatrix -> HeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsHeCRSMatrix<MA>::value
                 && IsHeMatrix<MB>::value,
         void>::Type
copy(const MA &A, MB &&B)
{
    typedef typename MA::IndexType    IndexType;
    typedef typename MA::ElementType  ElementType;

    if (B.dim()!=A.dim()) {
        ASSERT(B.dim()==0);
        B.resize(A.dim(), A.upLo(), A.indexBase(), ElementType(0));
        B.upLo() = A.upLo();
    }

    ASSERT(A.upLo()==B.upLo());

    const auto &rows = A.engine().rows();
    const auto &cols = A.engine().cols();
    const auto &vals = A.engine().values();

    for (IndexType i=rows.firstIndex(); i<rows.lastIndex(); ++i) {
        for (IndexType k=rows(i); k<rows(i+1); ++k) {
            B(i,cols(k)) = vals(k);
        }
    }
}

//-- copy: SyCRSMatrix -> SyMatrix
template <typename MA, typename MB>
typename RestrictTo<IsSyCRSMatrix<MA>::value
                 && IsSyMatrix<MB>::value,
         void>::Type
copy(const MA &A, MB &&B)
{
    typedef typename MA::IndexType    IndexType;
    typedef typename MA::ElementType  ElementType;

    if (B.dim()!=A.dim()) {
        ASSERT(B.dim()==0);
        B.resize(A.dim(), A.upLo(), A.indexBase(), ElementType(0));
        B.upLo() = A.upLo();
    }

    ASSERT(A.upLo()==B.upLo());

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
