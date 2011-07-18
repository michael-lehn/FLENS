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

#include <flens/storage/storageinfo.h>

namespace flens { namespace blas {

//-- product type: GeneralMatrix - Vector products -----------------------------
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
void
mv(cxxblas::Transpose trans,
   const ALPHA &alpha, const GeneralMatrix<MA> &A, const Vector<VX> &x,
   const BETA &beta, Vector<VY> &y)
{
    mv(trans, alpha, A.impl(), x.impl(), beta, y.impl());
}

//-- gemv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
void
mv(cxxblas::Transpose trans,
   const ALPHA &alpha, const GeMatrix<MA> &A, const DenseVector<VX> &x,
   const BETA &beta, DenseVector<VY> &y)
{
    FLENS_CLOSURELOG_ADD_ENTRY_GEMV(trans, alpha, A, x, beta, y);

    ASSERT(ADDRESS(y)!=ADDRESS(x));
    ASSERT(x.length()==((trans==cxxblas::NoTrans) ? A.numCols()
                                                  : A.numRows()));

    typedef typename GeMatrix<MA>::IndexType IndexType;
    IndexType yLength = (trans==cxxblas::NoTrans) ? A.numRows()
                                                  : A.numCols();

    ASSERT((beta==BETA(0)) || (y.length()==yLength));

    if (y.length()!=yLength) {
        y.engine().resize(yLength, 0);
    }

#   ifdef HAVE_CXXBLAS_GEMV
    cxxblas::gemv(StorageInfo<MA>::Order,
                  trans,
                  A.numRows(), A.numCols(),
                  alpha,
                  A.engine().data(), A.engine().leadingDimension(),
                  x.engine().data(), x.engine().stride(),
                  beta,
                  y.engine().data(), y.engine().stride());
#   else
    ASSERT(0);
#   endif

    FLENS_CLOSURELOG_END_ENTRY
}

//-- product type: HermitianMatrix - Vector products ---------------------------
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
void
mv(const ALPHA &alpha, const HermitianMatrix<MA> &A, const Vector<VX> &x,
   const BETA &beta, Vector<VY> &y)
{
    mv(alpha, A.impl(), x.impl(), beta, y.impl());
}

//-- hemv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
void
mv(const ALPHA &alpha, const HeMatrix<MA> &A, const DenseVector<VX> &x,
   const BETA &beta, DenseVector<VY> &y)
{
    ASSERT(ADDRESS(y)!=ADDRESS(x));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==BETA(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        y.engine().resize(A.dim(), 0);
    }

#   ifdef HAVE_CXXBLAS_HEMV
    cxxblas::hemv(StorageInfo<MA>::Order, A.upLo(),
                  A.dim(),
                  alpha,
                  A.engine().data(), A.engine().leadingDimension(),
                  x.engine().data(), x.engine().stride(),
                  beta,
                  y.engine().data(), y.engine().stride());
#   else
    ASSERT(0);
#   endif
}

//-- product type: SymmetricMatrix - Vector products ---------------------------
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
void
mv(const ALPHA &alpha, const SymmetricMatrix<MA> &A, const Vector<VX> &x,
   const BETA &beta, Vector<VY> &y)
{
    mv(alpha, A.impl(), x.impl(), beta, y.impl());
}

//-- symv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
void
mv(const ALPHA &alpha, const SyMatrix<MA> &A, const DenseVector<VX> &x,
   const BETA &beta, DenseVector<VY> &y)
{
    ASSERT(ADDRESS(y)!=ADDRESS(x));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==BETA(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        y.engine().resize(A.dim(), 0);
    }

#   ifdef HAVE_CXXBLAS_SYMV
    cxxblas::symv(StorageInfo<MA>::Order, A.upLo(),
                  A.dim(),
                  alpha,
                  A.engine().data(), A.engine().leadingDimension(),
                  x.engine().data(), x.engine().stride(),
                  beta,
                  y.engine().data(), y.engine().stride());
#   else
    ASSERT(0);
#   endif
}

//-- product type: TriangularMatrix - Vector products --------------------------
template <typename MA, typename VX>
void
mv(cxxblas::Transpose trans,  const TriangularMatrix<MA> &A, Vector<VX> &x)
{
    mv(trans, A.impl(), x.impl());
}

//-- trmv
template <typename MA, typename VX>
void
mv(cxxblas::Transpose trans, const TrMatrix<MA> &A, DenseVector<VX> &x)
{
    ASSERT(x.length()==A.dim());
#   ifdef HAVE_CXXBLAS_TRMV
    cxxblas::trmv(StorageInfo<MA>::Order, A.upLo(),
                  trans, A.diag(),
                  A.dim(),
                  A.engine().data(), A.engine().leadingDimension(),
                  x.engine().data(), x.engine().stride());
#   else
    ASSERT(0);
#   endif
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL2_MV_TCC
