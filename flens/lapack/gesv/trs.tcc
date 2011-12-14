/*
 *   Copyright (c) 2011, Michael Lehn
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

/* Besed on
 *
       SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_GESV_TRS_TCC
#define FLENS_LAPACK_GESV_TRS_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
template <typename MA, typename VP, typename MB>
void
trs_generic(Transpose trans, const GeMatrix<MA> &A, const DenseVector<VP> &piv,
            GeMatrix<MB> &B)
{
    using lapack::laswp;

    typedef typename GeMatrix<MA>::IndexType    IndexType;
    typedef typename GeMatrix<MA>::ElementType  T;

    const IndexType n       = A.numCols();
    const IndexType nRhs    = B.numCols();

//
//  Quick return if possible
//
    if ((n==0) || (nRhs==0)) {
        return;
    }

    if ((trans==NoTrans) || (trans==Conj)) {
//
//      Solve A * X = B.
//
//      Apply row interchanges to the right hand sides.
//
        laswp(B, piv);
//
//      Solve L*X = B, overwriting B with X.
//
        blas::sm(Left, trans, T(1), A.lowerUnit(), B);
//
//      Solve U*X = B, overwriting B with X.
//
        blas::sm(Left, trans, T(1), A.upper(), B);
    } else {
//
//      Solve A' * X = B.
//
//      Solve U'*X = B, overwriting B with X.
//
        blas::sm(Left, trans, T(1), A.upper(), B);
//
//      Solve L'*X = B, overwriting B with X.
//
        blas::sm(Left, trans, T(1), A.lowerUnit(), B);
//
//      Apply row interchanges to the solution vectors.
//
        laswp(B, piv.reverse());
    }
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename MA, typename VP, typename MB>
void
trs_native(Transpose trans, const GeMatrix<MA> &A, const DenseVector<VP> &piv,
           GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::ElementType ElementType;

    const char       TRANS = getF77LapackChar(trans);
    const INTEGER    N     = A.numRows();
    const INTEGER    NRHS  = B.numCols();
    const INTEGER    LDA   = A.leadingDimension();
    const INTEGER    LDB   = B.leadingDimension();
    INTEGER          INFO;


    if (IsSame<ElementType, double>::value) {
        LAPACK_IMPL(dgetrs)(&TRANS,
                            &N,
                            &NRHS,
                            A.data(),
                            &LDA,
                            piv.data(),
                            B.data(),
                            &LDB,
                            &INFO);
    } else {
        ASSERT(0);
    }
    ASSERT(INFO==0);
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename MA, typename VP, typename MB>
void
trs(Transpose trans, const GeMatrix<MA> &A, const DenseVector<VP> &piv,
    GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());

    const IndexType n = A.numRows();

    ASSERT(piv.firstIndex()==1);
    ASSERT(piv.length()==n);

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==n);
#   endif

//
//  Make copies of output arguments
//
    typename GeMatrix<MB>::NoView  B_org   = B;
//
//  Call implementation
//
    trs_generic(trans, A, piv, B);
//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    typename GeMatrix<MB>::NoView  B_generic   = B;

    B   = B_org;

    trs_native(trans, A, piv, B);

    bool failed = false;
    if (! isIdentical(B_generic, B, "B_generic", "B")) {
        std::cerr << "CXXLAPACK: B_generic = " << B_generic << std::endl;
        std::cerr << "F77LAPACK: B = " << B << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    } else {
        // std::cerr << "passed: (ge)trs.tcc" << std::endl;
    }

#   endif
}

template <typename MA, typename VP, typename VB>
void
trs(Transpose trans, const GeMatrix<MA> &A, const DenseVector<VP> &piv,
    DenseVector<VB> &b)
{
    typedef typename DenseVector<VB>::ElementType  ElementType;
    typedef typename DenseVector<VB>::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = GeMatrix<MA>::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b, n);

    return trs(trans, A, piv, B);
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VP, typename MB>
void
trs(Transpose trans, const MA &A, const VP &piv, MB &&B)
{
    trs(trans, A, piv, B);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GESV_TRS_TCC
