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

namespace flens { namespace blas {

//== product type: GeneralMatrix - GeneralMatrix products

//-- forwarding ----------------------------------------------------------------
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(cxxblas::Transpose transA, cxxblas::Transpose transB,
   const ALPHA &alpha,
   const MA &A, const MB &B,
   const BETA &beta,
   MC &&C)
{
    mm(transA, transB, alpha, A, B, beta, C);
}

//-- common interface ----------------------------------------------------------
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

    if (MC::order!=MA::order) {
        transA = cxxblas::Transpose(transA ^ cxxblas::Trans);
    }
    if (MC::order!=MB::order) {
        transB = cxxblas::Transpose(transB ^ cxxblas::Trans);
    }

    ASSERT((beta==BETA(0)) || (C.numRows()==m));
    ASSERT((beta==BETA(0)) || (C.numCols()==n));

    if ((C.numRows()!=m) || (C.numCols()!=n)) {
        C.resize(m, n);
    }

#   ifdef HAVE_CXXBLAS_GEMM
    cxxblas::gemm(MC::order,
                  transA, transB,
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
}


//== product type: HermitianMatrix - GeneralMatrix products

//-- forwarding ----------------------------------------------------------------
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(cxxblas::Side side, const ALPHA &alpha, const MA &A, const MB &B,
   const BETA &beta, MC &&C)
{
    mm(side, alpha, A, B, beta, C);
}

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
    ASSERT(MC::order==MB::Oorder);
    if (side==cxxblas::Left) {
        ASSERT(A.dim()==B.numRows());
    } else {
        ASSERT(B.numCols()==A.dim());
    }
#   endif

    cxxblas::StorageUpLo upLo = (MC::order==MA::order)
                              ? A.upLo()
                              : cxxblas::StorageUpLo(! A.upLo());

    typedef typename GeMatrix<MC>::IndexType IndexType;
    IndexType m = (side==cxxblas::Left) ? A.dim() : B.numRows();
    IndexType n = (side==cxxblas::Left) ? B.numCols() : A.dim();
 
    ASSERT((beta==static_cast<BETA>(0)) || (C.numRows()==m));
    ASSERT((beta==static_cast<BETA>(0)) || (C.numCols()==n));
 
    if ((C.numRows()!=m) || (C.numCols()!=n)) {
        C.resize(m,n);
    }

#   ifdef HAVE_CXXBLAS_HEMM
    cxxblas::hemm(MC::order, side,
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


//== product type: SymmetricMatrix - GeneralMatrix products

//-- forwarding ----------------------------------------------------------------
// -> is identical with forwarding of Hermitian Matrix - GeneralMatrix products

//-- common interface ----------------------------------------------------------
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
    ASSERT(MC::order==MB::order);
    if (side==cxxblas::Left) {
        ASSERT(A.dim()==B.numRows());
    } else {
        ASSERT(B.numCols()==A.dim());
    }
#   endif

    cxxblas::StorageUpLo upLo = (MC::order==MA::order)
                              ? A.upLo()
                              : cxxblas::StorageUpLo(! A.upLo());

    typedef typename GeMatrix<MC>::IndexType IndexType;
    IndexType m = (side==cxxblas::Left) ? A.dim() : B.numRows();
    IndexType n = (side==cxxblas::Left) ? B.numCols() : A.dim();

    ASSERT((beta==static_cast<BETA>(0)) || (C.numRows()==m));
    ASSERT((beta==static_cast<BETA>(0)) || (C.numCols()==n));

    if ((C.numRows()!=m) || (C.numCols()!=n)) {
        C.resize(m, n);
    }

#   ifdef HAVE_CXXBLAS_SYMM
    cxxblas::symm(MC::order, side,
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

//== product type: TriangularMatrix - GeneralMatrix products

//-- forwarding ----------------------------------------------------------------
template <typename ALPHA, typename MA, typename MB>
void
mm(cxxblas::Side side, cxxblas::Transpose transA,
   const ALPHA &alpha, const MA &A, MB &&B)
{
    mm(side, transA, alpha, A, B);
}

//-- common interface ----------------------------------------------------------
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
    ASSERT(MB::order==MA::order);
    if (side==cxxblas::Left) {
        assert(A.dim()==B.numRows());
    } else {
        assert(B.numCols()==A.dim());
    }
#   endif


#   ifdef HAVE_CXXBLAS_TRMM
    cxxblas::trmm(MB::order, side,
                  A.upLo(), transA, A.diag(),
                  B.numRows(), B.numCols(),
                  alpha,
                  A.data(), A.leadingDimension(),
                  B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL3_MM_TCC
