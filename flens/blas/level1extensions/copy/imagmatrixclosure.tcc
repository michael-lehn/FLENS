/*
 *   Copyright (c) 2012, Michael Lehn
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

#ifndef FLENS_BLAS_LEVEL1EXTENSIONS_COPY_IMAGMATRIXCLOSURE_TCC
#define FLENS_BLAS_LEVEL1EXTENSIONS_COPY_IMAGMATRIXCLOSURE_TCC 1

#include <flens/blas/level1extensions/copy/imagmatrixclosure.h>

namespace flens { namespace blas {

// ImagMatrixClosure -> GeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsGeMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const ImagMatrixClosure<MA> &A, MB &&B)
{
    const auto &Z = A._Z;

    ASSERT(trans==NoTrans || trans==Trans);

    if (B.numRows()==0 || B.numCols()==0) {
        if (trans==NoTrans) {
            B.resize(Z.numRows(), Z.numCols(), Z.firstRow(), Z.firstCol());
        } else {
            B.resize(Z.numCols(), Z.numRows(), Z.firstCol(), Z.firstRow());
        }
    }

#   ifndef NDEBUG
    if (trans==NoTrans) {
        ASSERT(Z.numRows()==B.numRows());
        ASSERT(Z.numCols()==B.numCols());
    } else {
        ASSERT(Z.numRows()==B.numCols());
        ASSERT(Z.numCols()==B.numRows());
    }
#   endif

    typedef typename RemoveRef<MB>::Type  MatrixB;
    typedef typename MatrixB::IndexType   IndexType;

    const Underscore<IndexType>  _;

    if (trans==NoTrans) {
        if (B.order()==RowMajor) {
            IndexType  i0 = Z.firstRow();
            IndexType  i1 = Z.lastRow();
            IndexType  I0 = B.firstRow();
            for (IndexType i=i0, I=I0; i<=i1; ++i, ++I) {
                B(I,_) = imag(Z(i,_));
            }
        } else {
            IndexType  j0 = Z.firstCol();
            IndexType  j1 = Z.lastCol();
            IndexType  J0 = B.firstCol();
            for (IndexType j=j0, J=J0; j<=j1; ++j, ++J) {
                B(_,J) = imag(Z(_,j));
            }
         }
    } else {
        if (B.order()==RowMajor) {
            IndexType  i0 = Z.firstRow();
            IndexType  i1 = Z.lastRow();
            IndexType  I0 = B.firstCol();
            for (IndexType i=i0, I=I0; i<=i1; ++i, ++I) {
                B(_,I) = imag(Z(i,_));
            }
        } else {
            IndexType  j0 = Z.firstCol();
            IndexType  j1 = Z.lastCol();
            IndexType  J0 = B.firstRow();
            for (IndexType j=j0, J=J0; j<=j1; ++j, ++J) {
                B(J,_) = imag(Z(_,j));
            }
        }
    }
}

// ImagConstMatrixClosure -> GeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsGeMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const ImagConstMatrixClosure<MA> &A, MB &&B)
{
    const auto &Z = A._Z;

    ASSERT(trans==NoTrans || trans==Trans);

    if (B.numRows()==0 || B.numCols()==0) {
        if (trans==NoTrans) {
            B.resize(Z.numRows(), Z.numCols(), Z.firstRow(), Z.firstCol());
        } else {
            B.resize(Z.numCols(), Z.numRows(), Z.firstCol(), Z.firstRow());
        }
    }

#   ifndef NDEBUG
    if (trans==NoTrans) {
        ASSERT(Z.numRows()==B.numRows());
        ASSERT(Z.numCols()==B.numCols());
    } else {
        ASSERT(Z.numRows()==B.numCols());
        ASSERT(Z.numCols()==B.numRows());
    }
#   endif

    typedef typename RemoveRef<MB>::Type  MatrixB;
    typedef typename MatrixB::IndexType   IndexType;

    const Underscore<IndexType>  _;

    if (trans==NoTrans) {
        if (B.order()==RowMajor) {
            IndexType  i0 = Z.firstRow();
            IndexType  i1 = Z.lastRow();
            IndexType  I0 = B.firstRow();
            for (IndexType i=i0, I=I0; i<=i1; ++i, ++I) {
                B(I,_) = imag(Z(i,_));
            }
        } else {
            IndexType  j0 = Z.firstCol();
            IndexType  j1 = Z.lastCol();
            IndexType  J0 = B.firstCol();
            for (IndexType j=j0, J=J0; j<=j1; ++j, ++J) {
                B(_,J) = imag(Z(_,j));
            }
         }
    } else {
        if (B.order()==RowMajor) {
            IndexType  i0 = Z.firstRow();
            IndexType  i1 = Z.lastRow();
            IndexType  I0 = B.firstCol();
            for (IndexType i=i0, I=I0; i<=i1; ++i, ++I) {
                B(_,I) = imag(Z(i,_));
            }
        } else {
            IndexType  j0 = Z.firstCol();
            IndexType  j1 = Z.lastCol();
            IndexType  J0 = B.firstRow();
            for (IndexType j=j0, J=J0; j<=j1; ++j, ++J) {
                B(J,_) = imag(Z(_,j));
            }
        }
    }
}

// GeMatrix -> ImagMatrixClosure
template <typename MA, typename MB>
void
copy(Transpose trans, const GeMatrix<MA> &A, ImagMatrixClosure<MB> &B)
{
    auto &Z = B._Z;

    ASSERT(trans==NoTrans || trans==Trans);

    if (Z.numRows()==0 || Z.numCols()==0) {
        if (trans==NoTrans) {
            Z.resize(A.numRows(), A.numCols(), A.firstRow(), A.firstCol());
        } else {
            Z.resize(A.numCols(), A.numRows(), A.firstCol(), A.firstRow());
        }
    }

#   ifndef NDEBUG
    if (trans==NoTrans) {
        ASSERT(Z.numRows()==A.numRows());
        ASSERT(Z.numCols()==A.numCols());
    } else {
        ASSERT(Z.numRows()==A.numCols());
        ASSERT(Z.numCols()==A.numRows());
    }
#   endif

    typedef typename MA::IndexType   IndexType;

    const Underscore<IndexType>  _;

    if (trans==NoTrans) {
        if (A.order()==RowMajor) {
            IndexType  i0 = Z.firstRow();
            IndexType  i1 = Z.lastRow();
            IndexType  I0 = A.firstRow();
            for (IndexType i=i0, I=I0; i<=i1; ++i, ++I) {
                imag(Z(i,_)) = A(I,_);
            }
        } else {
            IndexType  j0 = Z.firstCol();
            IndexType  j1 = Z.lastCol();
            IndexType  J0 = A.firstCol();
            for (IndexType j=j0, J=J0; j<=j1; ++j, ++J) {
                imag(Z(_,j)) = A(_,J);
            }
        }
    } else {
        if (A.order()==RowMajor) {
            IndexType  j0 = Z.firstCol();
            IndexType  j1 = Z.lastCol();
            IndexType  I0 = A.firstRow();
            for (IndexType j=j0, I=I0; j<=j1; ++j, ++I) {
                imag(Z(_,j)) = A(I,_);
            }
        } else {
            IndexType  i0 = Z.firstRow();
            IndexType  i1 = Z.lastRow();
            IndexType  J0 = A.firstCol();
            for (IndexType i=i0, J=J0; i<=i1; ++i, ++J) {
                imag(Z(i,_)) = A(_,J);
            }
        }
    }
}


} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1EXTENSIONS_COPY_IMAGMATRIXCLOSURE_TCC
