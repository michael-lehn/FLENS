/*
 *   Copyright (c) 2009, Michael Lehn
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

#ifndef FLENS_BLAS_EXTENSIONS_GEMATRIX_COPY_TCC
#define FLENS_BLAS_EXTENSIONS_GEMATRIX_COPY_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>

namespace flens { namespace blas {

//-- Densify Sparse Matrices ---------------------------------------------------

//-- copy: GeCCSMatrix -> GeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsGeCCSMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
copy(Transpose DEBUG_VAR(trans), const MA &A, MB &&B)
{
    ASSERT(trans==NoTrans);

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
copy(Transpose DEBUG_VAR(trans), const MA &A, MB &&B)
{
    ASSERT(trans==NoTrans);

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

//-- copy: GeCoordMatrix -> GeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsGeCoordMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
copy(Transpose DEBUG_VAR(trans), const MA &A, MB &&B)
{
    ASSERT(trans==NoTrans);

    typedef typename MA::ElementType  ElementType;

    B.resize(A.numRows(), A.numCols(),
             A.firstRow(), A.firstCol(),
             ElementType(0));

    const auto &coord = A.engine().coordVector();

    for (size_t k=0; k<coord.size(); ++k) {
        B(coord[k].row, coord[k].col) += coord[k].value;
    }
}

//-- Convenience Extensions ----------------------------------------------------

//-- copy: HeMatrix -> GeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsHeMatrix<MA>::value
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
        B.strictLower() = conjTrans(A.general().strictUpper());
    } else {
        B.lower() = A.general().lower();
        B.strictUpper() = conjTrans(A.general().strictLower());
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

} } // namespace blas, flens

#endif // FLENS_BLAS_EXTENSIONS_GEMATRIX_COPY_TCC
