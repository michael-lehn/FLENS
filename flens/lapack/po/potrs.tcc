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
       SUBROUTINE ZPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_PO_POTRS_TCC
#define FLENS_LAPACK_PO_POTRS_TCC 1

#include <cxxstd/algorithm.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- potrs [real variant] ------------------------------------------------------

template <typename MA, typename MB>
void
potrs_impl(const SyMatrix<MA> &A, GeMatrix<MB> &B)
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

//-- potrs [complex variant] ---------------------------------------------------

template <typename MA, typename MB>
void
potrs_impl(const HeMatrix<MA> &A, GeMatrix<MB> &B)
{
    using std::isnan;
    using std::sqrt;

    typedef typename HeMatrix<MA>::ElementType  ElementType;
    typedef typename HeMatrix<MA>::IndexType    IndexType;


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
//      Solve A*X = B where A = U**H *U.
//
//      Solve U**H *X = B, overwriting B with X.
//
        blas::sm(Left, ConjTrans, One, A.triangular(), B);
//
//      Solve U*X = B, overwriting B with X.
//
        blas::sm(Left, NoTrans, One, A.triangular(), B);
    } else {
//
//      Solve A*X = B where A = L*L**H.
//
//      Solve L*X = B, overwriting B with X.
//
        blas::sm(Left, NoTrans, One, A.triangular(), B);
//
//      Solve L**H *X = B, overwriting B with X.
//
        blas::sm(Left, ConjTrans, One, A.triangular(), B);
    }
}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- potrs [real variant] ------------------------------------------------------

template <typename MA, typename MB>
void
potrs_impl(const SyMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename SyMatrix<MA>::IndexType  IndexType;

    cxxlapack::potrs<IndexType>(getF77Char(A.upLo()),
                                A.dim(),
                                B.numCols(),
                                A.data(),
                                A.leadingDimension(),
                                B.data(),
                                B.leadingDimension());
}

//-- potrs [complex variant] ---------------------------------------------------

template <typename MA, typename MB>
void
potrs_impl(const HeMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename HeMatrix<MA>::IndexType  IndexType;

    cxxlapack::potrs<IndexType>(getF77Char(A.upLo()),
                                A.dim(),
                                B.numCols(),
                                A.data(),
                                A.leadingDimension(),
                                B.data(),
                                B.leadingDimension());
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- potrs [real/complex variant] ----------------------------------------------

template <typename MA, typename MB>
typename RestrictTo<(IsRealSyMatrix<MA>::value
                  && IsRealGeMatrix<MB>::value)
         ||         (IsHeMatrix<MA>::value
                  && IsComplexGeMatrix<MB>::value),
         void>::Type
potrs(const MA &A, MB &&B)
{
//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MB>::Type    MatrixB;
#   endif

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
    typename MatrixB::NoView B_org = B;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::potrs_impl(A, B);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of generic results
//
    typename MatrixB::NoView B_generic = B;
//
//  Restore output arguments
//
    B = B_org;

//
//  Compare results
//
    external::potrs_impl(A, B);

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

//-- potrs [variant if rhs is vector] ------------------------------------------

template <typename MA, typename VB>
typename RestrictTo<(IsSyMatrix<MA>::value || IsHeMatrix<MA>::value)
                 && IsDenseVector<VB>::value,
         void>::Type
potrs(const MA &A, VB &&b)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type     MatrixA;
    typedef typename RemoveRef<VB>::Type     VectorB;

//
//  Create matrix view from vector b and call above variant
//
    typedef typename VectorB::ElementType  ElementType;
    typedef typename VectorB::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = MatrixA::Engine::noViewOrder;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b);

    potrs(A, B);
}


} } // namespace lapack, flens

#endif // FLENS_LAPACK_PO_POTRS_TCC
