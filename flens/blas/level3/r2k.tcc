/*
 *   Copyright (c) 2010, Michael Lehn
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

#ifndef FLENS_BLAS_LEVEL3_R2K_TCC
#define FLENS_BLAS_LEVEL3_R2K_TCC

namespace flens { namespace blas {

//-- her2k
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
r2k(cxxblas::Transpose trans,
    const ALPHA &alpha,
    const GeMatrix<MA> &A, const GeMatrix<MB> &B,
    const BETA &beta,
    HeMatrix<MC> &C)
{
    typedef typename GeMatrix<MA>::IndexType IndexType;

    IndexType n = (trans==cxxblas::NoTrans) ? A.numRows()
                                            : A.numCols();
    IndexType k = (trans==cxxblas::NoTrans) ? A.numCols()
                                            : A.numRows();

    ASSERT(A.numRows()==B.numRows());
    ASSERT(A.numCols()==B.numCols());
    ASSERT((beta==static_cast<BETA>(0)) || (C.dim()==n));

    if (C.dim()!=n) {
        C.resize(n, n);
    }

#   ifdef HAVE_CXXBLAS_HER2K
    cxxblas::her2k(MC::order, C.upLo(), trans,
                   n, k,
                   alpha,
                   A.data(), A.leadingDimension(),
                   B.data(), B.leadingDimension(),
                   beta,
                   C.data(), C.leadingDimension());
#   else
    ASSERT(0);
#   endif
}

//-- syr2k
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
r2k(cxxblas::Transpose trans,
    const ALPHA &alpha,
    const GeMatrix<MA> &A, const GeMatrix<MB> &B,
    const BETA &beta,
    SyMatrix<MC> &C)
{
    typedef typename GeMatrix<MA>::IndexType IndexType;

    IndexType n = (trans==cxxblas::NoTrans) ? A.numRows()
                                            : A.numCols();
    IndexType k = (trans==cxxblas::NoTrans) ? A.numCols()
                                            : A.numRows();

    ASSERT(A.numRows()==B.numRows());
    ASSERT(A.numCols()==B.numCols());
    ASSERT((beta==static_cast<BETA>(0)) || (C.dim()==n));

    if (C.dim()!=n) {
        C.resize(n, n);
    }

#   ifdef HAVE_CXXBLAS_SYR2K
    cxxblas::syr2k(MC::order, C.upLo(), trans,
                   n, k,
                   alpha,
                   A.data(), A.leadingDimension(),
                   B.data(), B.leadingDimension(),
                   beta,
                   C.data(), C.leadingDimension());
#   else
    ASSERT(0);
#   endif
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL3_R2K_TCC
