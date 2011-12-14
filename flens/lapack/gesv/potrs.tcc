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
       SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_GESV_POTRS_TCC
#define FLENS_LAPACK_GESV_POTRS_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

#include <flens/lapack/interface/include/f77lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename MA, typename MB>
void
potrs_generic(const SyMatrix<MA> &A, GeMatrix<MB> &B)
{
    using std::isnan;
    using std::sqrt;

    typedef typename SyMatrix<MA>::ElementType  ElementType;
    typedef typename SyMatrix<MA>::IndexType    IndexType;


    const IndexType n    = A.dim();
    const IndexType nRhs = B.numCols();
    const bool upper     = (A.upLo()==Upper);

    const ElementType  One(1);
//
//  Quick return if possible
//
    if (n==0 || nRhs==0) {
        return;
    }
    if (upper) {
//
//      Solve A*X = B where A = U**T *U.
//
//      Solve U**T *X = B, overwriting B with X.
//
        blas::sm(Left, Trans, One, A.triangular(), B);
//
//      Solve U*X = B, overwriting B with X.
//
        blas::sm(Left, NoTrans, One, A.triangular(), B);
    } else {
//
//      Solve A*X = B where A = L*L**T.
//
//      Solve L*X = B, overwriting B with X.
//
        blas::sm(Left, NoTrans, One, A.triangular(), B);
//
//      Solve L**T *X = B, overwriting B with X.
//
        blas::sm(Left, Trans, One, A.triangular(), B);
    }
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename MA, typename MB>
void
potrs_native(const SyMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename SyMatrix<MA>::ElementType  T;

    const char       UPLO = char(A.upLo());
    const INTEGER    N    = A.dim();
    const INTEGER    NRHS = B.numCols();
    const INTEGER    LDA  = A.leadingDimension();
    const INTEGER    LDB  = B.leadingDimension();
    INTEGER          INFO;

    if (IsSame<T, DOUBLE>::value) {
        LAPACK_IMPL(dpotrs)(&UPLO, &N, &NRHS,
                            A.data(), &LDA,
                            B.data(), &LDB,
                            &INFO);
    } else {
        ASSERT(0);
    }
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename MA, typename MB>
void
potrs(const SyMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename SyMatrix<MA>::IndexType    IndexType;

//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);

    ASSERT(A.dim()==B.numRows());

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MB>::NoView B_org = B;
#   endif

//
//  Call implementation
//
    potrs_generic(A, B);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of generic results
//
    typename GeMatrix<MB>::NoView B_generic = B;
//
//  Restore output arguments
//
    B = B_org;

//
//  Compare results
//
    potrs_native(A, B);

    bool failed = false;
    if (! isIdentical(B_generic, B, "B_generic", "B")) {
        std::cerr << "CXXLAPACK: B_generic = " << B_generic << std::endl;
        std::cerr << "F77LAPACK: B = " << B << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename MB>
void
potrs(const MA &A, MB &&B)
{
    CHECKPOINT_ENTER;
    potrs(A, B);
    CHECKPOINT_LEAVE;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GESV_POTRS_TCC
