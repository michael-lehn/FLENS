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

/* Based on
 *
       SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 *
 *  -- LAPACK driver routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_GESV_SV_TCC
#define FLENS_LAPACK_GESV_SV_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
template <typename MA, typename VP, typename MB>
typename GeMatrix<MA>::IndexType
sv_generic(GeMatrix<MA> &A, DenseVector<VP> &piv, GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType IndexType;

    IndexType info = 0;

//
//  Compute the LU factorization of A.
//
    info = trf(A, piv);
    if (info==0) {
//
//      Solve the system A*X = B, overwriting B with X.
//
        trs(NoTrans, A, piv, B);
    }
    return info;
}

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename VP, typename MB>
typename GeMatrix<MA>::IndexType
sv(GeMatrix<MA> &A, DenseVector<VP> &piv, GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::gesv<IndexType>(A.numRows(),
                                                B.numCols(),
                                                A.data(),
                                                A.leadingDimension(),
                                                piv.data(),
                                                B.data(),
                                                B.leadingDimension());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename MA, typename VP, typename MB>
typename GeMatrix<MA>::IndexType
sv(GeMatrix<MA> &A, DenseVector<VP> &piv, GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType    IndexType;
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());
    ASSERT((piv.inc()>0 && piv.firstIndex()==1)
        || (piv.inc()<0 && piv.firstIndex()==A.numRows()));
    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==A.numRows());
#   endif
//
//  Make copies of output arguments
//
    typename GeMatrix<MA>::NoView    A_org   = A;
    typename DenseVector<VP>::NoView piv_org = piv;
    typename GeMatrix<MB>::NoView    B_org   = B;
//
//  Call implementation
//
    IndexType info = sv_generic(A, piv, B);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename GeMatrix<MA>::NoView    A_generic   = A;
    typename DenseVector<VP>::NoView piv_generic = piv;
    typename GeMatrix<MB>::NoView    B_generic   = B;

    A   = A_org;
    piv = piv_org;
    B   = B_org;

    IndexType _info = external::sv<IndexType>(A, piv, B);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "A_org = " << A_org << std::endl;
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(piv_generic, piv, "piv_generic", "piv")) {
        std::cerr << "piv_org = " << piv_org << std::endl;
        std::cerr << "CXXLAPACK: piv_generic = " << piv_generic << std::endl;
        std::cerr << "F77LAPACK: piv = " << piv << std::endl;
        failed = true;
    }

    if (! isIdentical(B_generic, B, "B_generic", "B")) {
        std::cerr << "B_org = " << B_org << std::endl;
        std::cerr << "CXXLAPACK: B_generic = " << B_generic << std::endl;
        std::cerr << "F77LAPACK: B = " << B << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, " info", "_info")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }
    if (failed) {
        ASSERT(0);
    }

#   endif

    return info;
}

template <typename MA, typename VP, typename VB>
typename GeMatrix<MA>::IndexType
sv(GeMatrix<MA> &A, DenseVector<VP> &piv, DenseVector<VB> &b)
{
    typedef typename DenseVector<VB>::ElementType  ElementType;
    typedef typename DenseVector<VB>::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = GeMatrix<MA>::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b, n);

    return sv(A, piv, B);
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VP, typename MB>
typename MA::IndexType
sv(MA &&A, VP &&piv, MB &&B)
{
    typedef typename MA::IndexType    IndexType;
    CHECKPOINT_ENTER;
    const IndexType info =  sv(A, piv, B);
    CHECKPOINT_LEAVE;
    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GESV_SV_TCC
