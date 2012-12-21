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
       SUBROUTINE DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
       SUBROUTINE ZPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 *
 *  -- LAPACK driver routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_PO_POSV_TCC
#define FLENS_LAPACK_PO_POSV_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- posv [real variant] -------------------------------------------------------

template <typename MA, typename MB>
typename SyMatrix<MA>::IndexType
posv_impl(SyMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MB>::IndexType  IndexType;

    IndexType info = potrf(A);

    if (info==0) {
        potrs(A, B);
    }
    return info;
}

//-- posv [complex variant] ----------------------------------------------------

template <typename MA, typename MB>
typename HeMatrix<MA>::IndexType
posv_impl(HeMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MB>::IndexType  IndexType;

    IndexType info = potrf(A);

    if (info==0) {
        potrs(A, B);
    }
    return info;
}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- posv [real variant] -------------------------------------------------------

template <typename MA, typename MB>
typename SyMatrix<MA>::IndexType
posv_impl(SyMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename SyMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::posv<IndexType>(getF77Char(A.upLo()),
                                                A.dim(),
                                                B.numCols(),
                                                A.data(),
                                                A.leadingDimension(),
                                                B.data(),
                                                B.leadingDimension());
    ASSERT(info>=0);
    return info;
}

//-- posv [complex variant] ----------------------------------------------------

template <typename MA, typename MB>
typename HeMatrix<MA>::IndexType
posv_impl(HeMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename HeMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::posv<IndexType>(getF77Char(A.upLo()),
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

//-- posv [real variant] -------------------------------------------------------

template <typename MA, typename MB>
typename RestrictTo<IsRealSyMatrix<MA>::value
                 && IsRealGeMatrix<MB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
posv(MA &&A, MB &&B)
{
    LAPACK_DEBUG_OUT("posv [real]");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename RemoveRef<MB>::Type    MatrixB;

//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);

    ASSERT(B.numRows()==A.dim());

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename MatrixA::NoView  A_org = A;
    typename MatrixB::NoView  B_org = B;
#   endif

//
//  Call implementation
//
    const IndexType info = LAPACK_SELECT::posv_impl(A, B);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename MatrixA::NoView  A_generic = A;
    typename MatrixB::NoView  B_generic = B;
    A = A_org;
    B = B_org;

    const IndexType _info = external::posv_impl(A, B);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "A_org = " << A_org << std::endl;
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
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


//
//  Complex variant
//
template <typename MA, typename MB>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsComplexGeMatrix<MB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
posv(MA &&A, MB &&B)
{
    LAPACK_DEBUG_OUT("posv [complex]");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename RemoveRef<MB>::Type    MatrixB;

//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);

    ASSERT(B.numRows()==A.dim());

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename MatrixA::NoView  A_org = A;
    typename MatrixB::NoView  B_org = B;
#   endif

//
//  Call implementation
//
    const IndexType info = LAPACK_SELECT::posv_impl(A, B);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename MatrixA::NoView  A_generic = A;
    typename MatrixB::NoView  B_generic = B;
    A = A_org;
    B = B_org;

    const IndexType _info = external::posv_impl(A, B);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "A_org = " << A_org << std::endl;
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
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



//-- posv [variant if rhs is vector] -----------------------------------------

template <typename MA, typename VB>
typename RestrictTo<(IsSyMatrix<MA>::value || IsHeMatrix<MA>::value)
                 && IsDenseVector<VB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
posv(MA &&A, VB &&b)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VB>::Type    VectorB;

    typedef typename VectorB::ElementType  ElementType;
    typedef typename VectorB::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = MatrixA::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b, n);

    return posv(A, B);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_PO_POSV_TCC
