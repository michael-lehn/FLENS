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
       SUBROUTINE DLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,
      $                   EQUED )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_EIG_LAQ_TCC
#define FLENS_LAPACK_EIG_LAQ_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
template <typename MA, typename VR, typename VC,
          typename ROWCOND, typename COLCOND,
          typename AMAX>
LAQ::Equilibration
laq_generic(GeMatrix<MA>            &A,
            const DenseVector<VR>   &r,
            const DenseVector<VC>   &c,
            const ROWCOND           &rowCond,
            const COLCOND           &colCond,
            const AMAX              &amax)
{
    using namespace LAQ;

    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const T One(1), Thresh(0.1);

    const Underscore<IndexType>  _;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();
//
//  Quick return if possible
//
    if (m==0 || n==0) {
        return None;
    }
//
//  Initialize LARGE and SMALL.
//
    const T small = lamch<T>(SafeMin) / lamch<T>(Precision);
    const T large = One / small;

    Equilibration  equed;

    if (rowCond>=Thresh && amax>=small && amax<=large) {
//
//  No row scaling
//
        if (colCond>=Thresh) {
//
//          No column scaling
//
            equed = None;
        } else {
//
//          Column scaling
//
            for (IndexType j=1; j<=n; ++j) {
                A(_,j) *= c(j);
            }
            equed = Column;
        }
    } else if (colCond>=Thresh) {
//
//      Row scaling, no column scaling
//
        for (IndexType j=1; j<=n; ++j) {
            for (IndexType i=1; i<=m; ++i) {
                A(i,j) *= r(i);
            }
        }
        equed = Row;
    } else {
//
//      Row and column scaling
//
        for (IndexType j=1; j<=n; ++j) {
            const T cj = c(j);
            for (IndexType i=1; i<=m; ++i) {
                A(i,j) *= cj*r(i);
            }
        }
        equed = Both;
    }
    return equed;
}


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename VR, typename VC,
          typename ROWCOND, typename COLCOND,
          typename AMAX>
LAQ::Equilibration
laq(GeMatrix<MA>            &A,
    const DenseVector<VR>   &r,
    const DenseVector<VC>   &c,
    const ROWCOND           &rowCond,
    const COLCOND           &colCond,
    const AMAX              &amax)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;
    char equed;

    cxxlapack::laqge<IndexType>(A.numRows(),
                                A.numCols(),
                                A.data(),
                                A.leadingDimension(),
                                r.data(),
                                c.data(),
                                rowCond,
                                colCond,
                                amax,
                                equed);

    return LAQ::Equilibration(equed);
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
template <typename MA, typename VR, typename VC,
          typename ROWCOND, typename COLCOND,
          typename AMAX>
LAQ::Equilibration
laq(GeMatrix<MA>            &A,
    const DenseVector<VR>   &r,
    const DenseVector<VC>   &c,
    const ROWCOND           &rowCond,
    const COLCOND           &colCond,
    const AMAX              &amax)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    ASSERT(r.firstIndex()==1);
    ASSERT(r.length()==m);

    ASSERT(c.firstIndex()==1);
    ASSERT(c.length()==n);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MA>::NoView  A_org  = A;
#   endif

//
//  Call implementation
//
    const auto equed = laq_generic(A, r, c, rowCond, colCond, amax);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename GeMatrix<MA>::NoView A_generic = A;
//
//  restore output arguments
//
    A = A_org;

//
//  Compare generic results with results from the native implementation
//
    const auto equed_ = external::laq(A, r, c, rowCond, colCond, amax);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = "
                  << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (equed!=equed_) {
        std::cerr << "CXXLAPACK: equed = " << equed << std::endl;
        std::cerr << "F77LAPACK: equed_ = " << equed_ << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    } else {
//        std::cerr << "passed: laq.tcc" << std::endl;
    }
#   endif

    return equed;
}


//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VR, typename VC,
          typename ROWCOND, typename COLCOND,
          typename AMAX>
LAQ::Equilibration
laq(MA             &&A,
    const VR       &r,
    const VC       &c,
    const ROWCOND  &rowCond,
    const COLCOND  &colCond,
    const AMAX     &amax)
{
    CHECKPOINT_ENTER;
    const auto equed = laq(A, r, c, rowCond, colCond, amax);
    CHECKPOINT_LEAVE;

    return equed;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_LAQ_TCC
