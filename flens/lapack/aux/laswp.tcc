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
 *   Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.
 *
 *   $COPYRIGHT$
 *
 *   Additional copyrights may follow
 *
 *   $HEADER$
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions are
 *   met:
 *
 *   - Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer. 
 *
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer listed
 *     in this license in the documentation and/or other materials
 *     provided with the distribution.
 *
 *   - Neither the name of the copyright holders nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
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

#ifndef FLENS_LAPACK_AUX_LASWP_TCC
#define FLENS_LAPACK_AUX_LASWP_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename MA, typename VP>
void
laswp_generic(GeMatrix<MA> &A, const DenseVector<VP> &piv)
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

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename MA, typename VP>
void
laswp_native(GeMatrix<MA> &_A, const DenseVector<VP> &piv)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MA>::ElementType  T;

    const INTEGER    N = _A.numCols();
    T               *A = _A.data();
    const INTEGER    LDA = _A.leadingDimension();
    const INTEGER    K1 = min(piv.firstIndex(), piv.lastIndex());
    const INTEGER    K2 = max(piv.firstIndex(), piv.lastIndex());

    // set pointer IPIV such that IPIV[K1] points to piv.data()
    // (assumes an index base of 1)
    const INTEGER   *IPIV = piv.data() - K1 + 1;
    INTEGER          INCX = piv.inc();

    if (IsSame<T, DOUBLE>::value) {
        LAPACK_IMPL(dlaswp)(&N, A, &LDA, &K1, &K2, IPIV, &INCX);
    } else {
        ASSERT(0);
    }
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename MA, typename VP>
void
laswp(GeMatrix<MA> &A, const DenseVector<VP> &piv)
{
    LAPACK_DEBUG_OUT("laswp");

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
    laswp_generic(A, piv);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    laswp_native(_A, piv);

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
    laswp(A, piv);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_LASWP_TCC
