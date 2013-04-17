/*
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

/* Based on
 *
       SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
       SUBROUTINE ZTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 *
 */

#ifndef FLENS_LAPACK_TR_TRS_TCC
#define FLENS_LAPACK_TR_TRS_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- (tr)trs [real and complex variant] ----------------------------------------

template <typename MA, typename MB>
typename TrMatrix<MA>::IndexType
trs_impl(Transpose trans, const TrMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename TrMatrix<MA>::IndexType    IndexType;
    typedef typename TrMatrix<MA>::ElementType  T;

    const IndexType n       = A.dim();

    const T  Zero(0), One(1);

    IndexType info = 0;
//
//  Quick return if possible
//
    if (n==0) {
        return info;
    }
//
//  Check for singularity.
//
    if (A.diag()!=Unit) {
        for (info=1; info<=n; ++info) {
            if (A(info,info)==Zero) {
                return info;
            }
        }
    }
    info = 0;
//
//  Solve A * x = b  or  A**T * x = b.
//
    blas::sm(Left, trans, One, A, B);

    return info;
}

} // namespace generic


//== interface for external lapack =============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (tr)trs [real and complex variant] ----------------------------------------

template <typename MA, typename MB>
typename TrMatrix<MA>::IndexType
trs_impl(Transpose trans, const TrMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename TrMatrix<MA>::IndexType  IndexType;

    IndexType info;
    info = cxxlapack::trtrs<IndexType>(getF77Char(A.upLo()),
                                       getF77Char(trans),
                                       getF77Char(A.diag()),
                                       A.dim(),
                                       B.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       B.data(),
                                       B.leadingDimension());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- (tr)trs [real and complex variant] ----------------------------------------

template <typename MA, typename MB>
typename RestrictTo<IsTrMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
trs(Transpose trans, const MA &A, MB &&B)
{
    LAPACK_DEBUG_OUT("(tr)trs [real/complex]");
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    const IndexType n = A.dim();

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==n);
#   endif

#   ifdef CHECK_CXXLAPACK

    typedef typename RemoveRef<MB>::Type    MatrixB;
//
//  Make copies of output arguments
//
    typename MatrixB::NoView  B_org   = B;
#   endif
//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::trs_impl(trans, A, B);
//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixB::NoView  B_generic   = B;

    B   = B_org;

    IndexType _info = external::trs_impl(trans, A, B);

    bool failed = false;
    if (! isIdentical(B_generic, B, "B_generic", "B")) {
        std::cerr << "CXXLAPACK: B_generic = " << B_generic << std::endl;
        std::cerr << "F77LAPACK: B = " << B << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, "info", "_info")) {
        std::cerr << "CXXLAPACK: info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    } else {
        // std::cerr << "passed: (tr)trs.tcc" << std::endl;
    }
#   endif

    return info;
}

//-- (tr)trs [variant if rhs is vector] ----------------------------------------

template <typename MA, typename VB>
typename RestrictTo<IsTrMatrix<MA>::value
                 && IsDenseVector<VB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
trs(Transpose trans, const MA &A, VB &&b)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VB>::Type    VectorB;

//
//  Create matrix view from vector b and call above variant
//
    typedef typename VectorB::ElementType  ElementType;
    typedef typename VectorB::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = MatrixA::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b, n);

    return trs(trans, A, B);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_TR_TRS_TCC
