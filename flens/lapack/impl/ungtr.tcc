/*
 *   Copyright (c) 2014, Michael Lehn
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
       SUBROUTINE ZUNGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_UNGTR_TCC
#define FLENS_LAPACK_IMPL_UNGTR_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- ungtr ---------------------------------------------------------------------

template <typename MA, typename VTAU, typename VWORK>
void
ungtr_impl(HeMatrix<MA>              &A,
           const DenseVector<VTAU>   &tau,
           DenseVector<VWORK>        &work)
{
    using std::max;

    typedef typename HeMatrix<MA>::ElementType       T;
    typedef typename HeMatrix<MA>::IndexType         IndexType;

    const T  Zero(0), One(1);

    const Underscore<IndexType> _;
    const IndexType n = A.numCols();

    auto Q = A.general();

    IndexType nb = (A.upLo()==Upper)
                 ? ilaenv<T>(1, "UNGQL", "", n-1, n-1, n-1)
                 : ilaenv<T>(1, "UNGQR", "", n-1, n-1, n-1);
    IndexType lWorkOpt = max(1,n-1)*nb;

//
//  Quick return if possible
//
    if (n==0) {
        work(1) = 1;
        return;
    }

    if (A.upLo()==Upper) {
//
//      Q was determined by a call to ZHETRD with UPLO = 'U'
//
//      Shift the vectors which define the elementary reflectors one
//      column to the left, and set the last row and column of Q to
//      those of the unit matrix
//
        for (IndexType j=1; j<=n-1; ++j) {
            for (IndexType i=1; i<=j-1; ++i) {
                Q(i,j) = Q(i,j+1);
            }
            Q(n,j) = Zero;
        }
        for (IndexType i=1; i<=n-1; ++i) {
            Q(i,n) = Zero;
        }
        Q(n,n) = One;
//
//      Generate Q(1:n-1,1:n-1)
//
        ungql(Q(_(1,n-1),_(1,n-1)), tau(_(1,n-1)), work);

    } else {
//
//      Q was determined by a call to ZHETRD with UPLO = 'L'.
//
//      Shift the vectors which define the elementary reflectors one
//      column to the right, and set the first row and column of Q to
//      those of the unit matrix
//
        for (IndexType j=n; j>=2; --j) {
            Q(1,j) = Zero;
            for (IndexType i=j+1; i<=n; ++i) {
                Q(i,j) = Q(i,j-1);
            }
        }
        Q(1,1) = One;
        for (IndexType i=2; i<=n; ++i) {
            Q(i,1) = Zero;
        }
        if (n>1) {
//
//          Generate Q(2:n,2:n)
//
            ungqr(Q(_(2,n),_(2,n)), tau(_(1,n-1)), work);
        }
    }
    work(1) = lWorkOpt;
}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- ungtr ---------------------------------------------------------------------

template <typename MA, typename VTAU, typename VWORK>
void
ungtr_impl(HeMatrix<MA>              &A,
           const DenseVector<VTAU>   &tau,
           DenseVector<VWORK>        &work)
{
    typedef typename HeMatrix<MA>::ElementType   ElementType;
    typedef typename HeMatrix<MA>::IndexType     IndexType;

    if (work.length()==0) {
        ElementType  WORK;
        IndexType    LWORK = -1;

        cxxlapack::ungtr<IndexType>(getF77Char(A.upLo()),
                                    A.dim(),
                                    A.data(),
                                    A.leadingDimension(),
                                    tau.data(),
                                    &WORK,
                                    LWORK);
        work.resize(IndexType(real(WORK)));
    }

    cxxlapack::ungtr<IndexType>(getF77Char(A.upLo()),
                                A.dim(),
                                    A.data(),
                                A.leadingDimension(),
                                tau.data(),
                                work.data(),
                                work.length());
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- ungtr ---------------------------------------------------------------------

template <typename MA, typename VTAU, typename VWORK>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsComplexDenseVector<VTAU>::value
                 && IsComplexDenseVector<VWORK>::value,
         void>::Type
ungtr(MA &&A, const VTAU &tau, VWORK &&work)
{
//
//  Remove references from rvalue types
//
#   if !defined(NDEBUG) || defined(CHECK_CXXLAPACK)
    typedef typename RemoveRef<MA>::Type    MatrixA;
#   endif

#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<VWORK>::Type VectorWork;
#   endif

//
//  Test the input parameters
//
#   ifndef NDEBUG
    typedef typename MatrixA::IndexType     IndexType;

    ASSERT(A.firstRow()==IndexType(1));
    ASSERT(A.firstCol()==IndexType(1));
    ASSERT(tau.firstIndex()==IndexType(1));
    ASSERT(tau.length()==A.numRows()-1);
    ASSERT((work.length()==0) || (work.length()>=tau.length()));
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView        A_org      = A;
    typename VectorWork::NoView     work_org   = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::ungtr_impl(A, tau, work);

#   ifdef CHECK_CXXLAPACK
//
//  Restore output arguments
//
    typename MatrixA::NoView        A_generic      = A;
    typename VectorWork::NoView     work_generic   = work;

    A    = A_org;
    work = work_org;
//
//  Compare results
//
    external::ungtr_impl(A, tau, work);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "char(A.upLo()) = " << char(A.upLo()) << std::endl;
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: ungtr.tcc" << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: ungtr.tcc" << std::endl;
    }
#   endif
}

//-- ungtr [Variant with temporary workspace] ----------------------------------

template <typename MA, typename VTAU>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsComplexDenseVector<VTAU>::value,
         void>::Type
ungtr(MA &&A, const VTAU &tau)
{
    typedef typename RemoveRef<MA>::Type::Vector  WorkVector;

    WorkVector  work;
    ungtr(A, tau, work);
}


} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_UNGTR_TCC
