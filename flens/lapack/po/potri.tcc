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
       SUBROUTINE DPOTRI( UPLO, N, A, LDA, INFO )
       SUBROUTINE ZPOTRI( UPLO, N, A, LDA, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_PO_POTRI_TCC
#define FLENS_LAPACK_PO_POTRI_TCC 1

#include <cxxstd/algorithm.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- potri [real variant] ------------------------------------------------------

template <typename MA>
typename SyMatrix<MA>::IndexType
potri_impl(SyMatrix<MA> &A)
{
    typedef typename SyMatrix<MA>::IndexType    IndexType;
    const IndexType n = A.dim();
    IndexType info = 0;
//
//  Quick return if possible
//
    if (n==0) {
        return 0;
    }
//
//  Invert the triangular Cholesky factor U or L.
//
    auto T = A.triangular();

    info = tri(T);
    if (info==0) {
//
//      Form inv(U) * inv(U)**T or inv(L)**T * inv(L).
//
        lauum(T);
    }
    return info;
}

//-- potri [complex variant] ---------------------------------------------------

template <typename MA>
typename HeMatrix<MA>::IndexType
potri_impl(HeMatrix<MA> &A)
{
    typedef typename HeMatrix<MA>::IndexType    IndexType;
    const IndexType n = A.dim();
    IndexType info = 0;
//
//  Quick return if possible
//
    if (n==0) {
        return 0;
    }
//
//  Invert the triangular Cholesky factor U or L.
//
    auto T = A.triangular();

    info = tri(T);
    if (info==0) {
//
//      Form inv(U) * inv(U)**T or inv(L)**T * inv(L).
//
        lauum(T);
    }
    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- potri [real variant] ------------------------------------------------------

template <typename MA>
typename SyMatrix<MA>::IndexType
potri_impl(SyMatrix<MA> &A)
{
    typedef typename SyMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::potri<IndexType>(getF77Char(A.upLo()),
                                                 A.dim(),
                                                 A.data(),
                                                 A.leadingDimension());
    ASSERT(info>=0);
    return info;
}

//-- potri [complex variant] ---------------------------------------------------

template <typename MA>
typename HeMatrix<MA>::IndexType
potri_impl(HeMatrix<MA> &A)
{
    typedef typename HeMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::potri<IndexType>(getF77Char(A.upLo()),
                                                 A.dim(),
                                                 A.data(),
                                                 A.leadingDimension());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//-- potri [real variant] ------------------------------------------------------

template <typename MA>
typename RestrictTo<IsRealSyMatrix<MA>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
potri(MA &&A)
{
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
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename MatrixA::NoView  A_org = A;
#   endif

//
//  Call implementation
//
    const IndexType info = LAPACK_SELECT::potri_impl(A);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of generic results
//
    typename MatrixA::NoView  A_generic = A;
//
//  restore output parameters
//
    A = A_org;
//
//  Compare results
//
    const IndexType _info = external::potri_impl(A);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "A_org =" << A_org << std::endl;
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
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

//-- potri [complex variant] ---------------------------------------------------

template <typename MA>
typename RestrictTo<IsHeMatrix<MA>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
potri(MA &&A)
{
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
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename MatrixA::NoView  A_org = A;
#   endif

//
//  Call implementation
//
    const IndexType info = LAPACK_SELECT::potri_impl(A);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of generic results
//
    typename MatrixA::NoView  A_generic = A;
//
//  restore output parameters
//
    A = A_org;
//
//  Compare results
//
    const IndexType _info = external::potri_impl(A);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "A_org =" << A_org << std::endl;
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
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

} } // namespace lapack, flens

#endif // FLENS_LAPACK_PO_POTRI_TCC
