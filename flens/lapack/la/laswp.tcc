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

/*
       SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LASWP_TCC
#define FLENS_LAPACK_LA_LASWP_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename MA, typename VP>
void
laswp_impl(GeMatrix<MA> &A, const DenseVector<VP> &piv)
{
    typedef typename GeMatrix<MA>::IndexType IndexType;
    typedef Range<IndexType>                 Range;

    const IndexType n = A.numCols();

//
//  Interchange row i with row piv(i) for each of rows K1 through K2.
//
    const IndexType pBeg = piv.firstIndex();
    const IndexType pEnd = piv.endIndex();
    const IndexType pInc = piv.inc();

    const IndexType bs = 32;
    const IndexType nbs = (n/bs)*bs;

    if (nbs!=0) {
        for (IndexType j=1; j<=nbs; j+=bs) {
            const Range cols(j,j+bs-1);
            for (IndexType i=pBeg; i!=pEnd; i+=pInc) {
                const IndexType iP = piv(i);
                if (iP!=i) {
                    blas::swap(A(i, cols), A(iP, cols));
                }
            }
        }
    }
    if (nbs!=n) {
        const Range cols(nbs+1,n);
        for (IndexType i=pBeg; i!=pEnd; i+=pInc) {
            const IndexType iP = piv(i);
            if (iP!=i) {
                blas::swap(A(i, cols), A(iP, cols));
            }
        }
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename VP>
void
laswp_impl(GeMatrix<MA> &A, const DenseVector<VP> &piv)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MA>::IndexType  IndexType;

    const IndexType k1 = min(piv.firstIndex(), piv.lastIndex());
    const IndexType k2 = max(piv.firstIndex(), piv.lastIndex());

    cxxlapack::laswp<IndexType>(A.numCols(), A.data(), A.leadingDimension(),
                                k1, k2,
                                piv.data()-k1+1, piv.inc());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename MA, typename VP>
void
laswp(GeMatrix<MA> &A, const DenseVector<VP> &piv)
{
//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT((piv.inc()>0 && piv.firstIndex()>=A.firstRow())
        || (piv.inc()<0 && piv.firstIndex()<=A.numRows()));
    ASSERT((piv.inc()>0 && piv.lastIndex()<=A.numRows())
        || (piv.inc()<0 && piv.lastIndex()>=A.firstRow()));

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MA>::NoView       org_A   = A;
    typename GeMatrix<MA>::NoView       _A      = A;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::laswp_impl(A, piv);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    external::laswp_impl(_A, piv);

    bool failed = false;
    if (! isIdentical(A, _A, " A", "_A")) {
        std::cerr << "CXXLAPACK:  A = " << A << std::endl;
        std::cerr << "F77LAPACK: _A = " << _A << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "piv.firstIndex() = " << piv.firstIndex() << std::endl;
        std::cerr << "piv.lastIndex() = " << piv.lastIndex() << std::endl;
        std::cerr << "piv.inc() = " << piv.inc() << std::endl;
        std::cerr << "piv = " << piv << std::endl;

        std::cerr << "org_A = " << org_A << std::endl;
        ASSERT(0);
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VP>
void
laswp(MA &&A, const VP &&piv)
{
    CHECKPOINT_ENTER;
    laswp(A, piv);
    CHECKPOINT_LEAVE;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LASWP_TCC
