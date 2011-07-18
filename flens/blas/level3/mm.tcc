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

#include <flens/storage/storageinfo.h>

namespace flens { namespace blas {

//-- product type: GeneralMatrix - GeneralMatrix products ----------------------
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(cxxblas::Transpose transA, cxxblas::Transpose transB,
   const ALPHA &alpha,
   const GeneralMatrix<MA> &A, const GeneralMatrix<MB> &B,
   const BETA &beta,
   GeneralMatrix<MC> &C)
{
    mm(transA, transB, alpha, A.impl(), B.impl(), beta, C.impl());
}

//-- gemm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(cxxblas::Transpose transA, cxxblas::Transpose transB,
   const ALPHA &alpha,
   const GeMatrix<MA> &A, const GeMatrix<MB> &B,
   const BETA &beta,
   GeMatrix<MC> &C)
{
#   ifndef NDEBUG
    int kA = (transA==cxxblas::NoTrans) ? A.numCols() : A.numRows();
    int kB = (transB==cxxblas::NoTrans) ? B.numRows() : B.numCols();
    ASSERT(kA==kB);
#   endif

    typedef typename GeMatrix<MC>::IndexType IndexType;
    IndexType m = (transA==cxxblas::NoTrans) ? A.numRows() : A.numCols();
    IndexType n = (transB==cxxblas::NoTrans) ? B.numCols() : B.numRows();
    IndexType k = (transA==cxxblas::NoTrans) ? A.numCols() : A.numRows();

    if (StorageInfo<MC>::Order!=StorageInfo<MA>::Order) {
        transA = cxxblas::Transpose(transA ^ cxxblas::Trans);
    }
    if (StorageInfo<MC>::Order!=StorageInfo<MB>::Order) {
        transB = cxxblas::Transpose(transB ^ cxxblas::Trans);
    }

    ASSERT((beta==BETA(0)) || (C.numRows()==m));
    ASSERT((beta==BETA(0)) || (C.numCols()==n));

    if ((C.numRows()!=m) || (C.numCols()!=n)) {
        C.engine().resize(m, n);
    }

#   ifdef HAVE_CXXBLAS_GEMM
    cxxblas::gemm(StorageInfo<MC>::Order,
                  transA, transB,
                  C.numRows(),
                  C.numCols(),
                  k,
                  alpha,
                  A.engine().data(), A.engine().leadingDimension(),
                  B.engine().data(), B.engine().leadingDimension(),
                  beta,
                  C.engine().data(), C.engine().leadingDimension());
#   else
    ASSERT(0);
#   endif
}

//-- product type: HermitianMatrix - GeneralMatrix products --------------------
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(cxxblas::Side side,
   const ALPHA &alpha,
   const HermitianMatrix<MA> &A, const GeneralMatrix<MB> &B,
   const BETA &beta, GeneralMatrix<MC> &C)
{
    mm(side, alpha, A.impl(), B.impl(), beta, C.impl());
}

//-- hemm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(cxxblas::Side side,
   const ALPHA &alpha, const HeMatrix<MA> &A, const GeMatrix<MB> &B,
   const BETA &beta, GeMatrix<MC> &C)
{
#   ifndef NDEBUG
    ASSERT(StorageInfo<MC>::Order==StorageInfo<MB>::Order);
    if (side==cxxblas::Left) {
        ASSERT(A.dim()==B.numRows());
    } else {
        ASSERT(B.numCols()==A.dim());
    }
#   endif

    cxxblas::StorageUpLo upLo = (StorageInfo<MC>::Order==StorageInfo<MA>::Order)
                              ? A.upLo()
                              : cxxblas::StorageUpLo(! A.upLo());

    typedef typename GeMatrix<MC>::IndexType IndexType;
    IndexType m = (side==cxxblas::Left) ? A.dim() : B.numRows();
    IndexType n = (side==cxxblas::Left) ? B.numCols() : A.dim();
 
    ASSERT((beta==static_cast<BETA>(0)) || (C.numRows()==m));
    ASSERT((beta==static_cast<BETA>(0)) || (C.numCols()==n));
 
    if ((C.numRows()!=m) || (C.numCols()!=n)) {
        C.engine().resize(m,n);
    }

#   ifdef HAVE_CXXBLAS_HEMM
    cxxblas::hemm(StorageInfo<MC>::Order, side,
                  upLo,
                  C.numRows(), C.numCols(),
                  alpha,
                  A.engine().data(), A.engine().leadingDimension(),
                  B.engine().data(), B.engine().leadingDimension(),
                  beta,
                  C.engine().data(), C.engine().leadingDimension());
#    else
     ASSERT(0);
#    endif
}

//-- product type: SymmetricMatrix - GeneralMatrix products --------------------
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(cxxblas::Side side,
   const ALPHA &alpha,
   const SymmetricMatrix<MA> &A, const GeneralMatrix<MB> &B,
   const BETA &beta, GeneralMatrix<MC> &C)
{
    mm(side, alpha, A.impl(), B.impl(), beta, C.impl());
}

//-- symm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(cxxblas::Side side,
   const ALPHA &alpha, const SyMatrix<MA> &A, const GeMatrix<MB> &B,
   const BETA &beta, GeMatrix<MC> &C)
{
#   ifndef NDEBUG
    ASSERT(StorageInfo<MC>::Order==StorageInfo<MB>::Order);
    if (side==cxxblas::Left) {
        ASSERT(A.dim()==B.numRows());
    } else {
        ASSERT(B.numCols()==A.dim());
    }
#   endif

    cxxblas::StorageUpLo upLo = (StorageInfo<MC>::Order==StorageInfo<MA>::Order)
                              ? A.upLo()
                              : cxxblas::StorageUpLo(! A.upLo());

    typedef typename GeMatrix<MC>::IndexType IndexType;
    IndexType m = (side==cxxblas::Left) ? A.dim() : B.numRows();
    IndexType n = (side==cxxblas::Left) ? B.numCols() : A.dim();

    ASSERT((beta==static_cast<BETA>(0)) || (C.numRows()==m));
    ASSERT((beta==static_cast<BETA>(0)) || (C.numCols()==n));

    if ((C.numRows()!=m) || (C.numCols()!=n)) {
        C.engine().resize(m, n);
    }

#   ifdef HAVE_CXXBLAS_SYMM
    cxxblas::symm(StorageInfo<MC>::Order, side,
                  upLo,
                  C.numRows(), C.numCols(),
                  alpha,
                  A.engine().data(), A.engine().leadingDimension(),
                  B.engine().data(), B.engine().leadingDimension(),
                  beta,
                  C.engine().data(), C.engine().leadingDimension());
#    else
     ASSERT(0);
#    endif
}

//-- product type: TriangularMatrix - GeneralMatrix products -------------------
template <typename ALPHA, typename MA, typename MB>
void
mm(cxxblas::Side side,
   cxxblas::Transpose transA, const ALPHA &alpha,
   const TriangularMatrix<MA> &A,
   GeneralMatrix<MB> &B)
{
    mm(side, transA, alpha, A.impl(), B.impl());
}

//-- trmm
template <typename ALPHA, typename MA, typename MB>
void
mm(cxxblas::Side side,
   cxxblas::Transpose transA, const ALPHA &alpha, const TrMatrix<MA> &A,
   GeMatrix<MB> &B)
{
#   ifndef NDEBUG
    ASSERT(StorageInfo<MB>::Order==StorageInfo<MA>::Order);
    if (side==cxxblas::Left) {
        assert(A.dim()==B.numRows());
    } else {
        assert(B.numCols()==A.dim());
    }
#   endif


#   ifdef HAVE_CXXBLAS_TRMM
    cxxblas::trmm(StorageInfo<MB>::Order, side,
                  A.upLo(), transA, A.diag(),
                  B.numRows(), B.numCols(),
                  alpha,
                  A.engine().data(), A.engine().leadingDimension(),
                  B.engine().data(), B.engine().leadingDimension());
#   else
    ASSERT(0);
#   endif
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL3_MM_TCC
