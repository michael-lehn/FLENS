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

#ifndef FLENS_BLAS_EXTENSIONS_SYMATRIX_COPY_TCC
#define FLENS_BLAS_EXTENSIONS_SYMATRIX_COPY_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/level1/level1.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>

namespace flens { namespace blas {

//-- copy: SyMatrix -> SyMatrix
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

    if (A.upLo()==B.upLo()) {
        copy(NoTrans, A.triangular(), B.triangular());
    } else {
        copy(ConjTrans, A.triangular(), B.triangular());
    }
}

//-- Densify Sparse Matrices ---------------------------------------------------

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

#endif // FLENS_BLAS_EXTENSIONS_SYMATRIX_COPY_TCC
