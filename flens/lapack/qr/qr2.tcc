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

#ifndef FLENS_LAPACK_QR_QR2_TCC
#define FLENS_LAPACK_QR_QR2_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename MA, typename VTAU, typename VWORK>
void
qr2_generic(GeMatrix<MA> &A, DenseVector<VTAU> &tau, DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MA>::IndexType    IndexType;
    typedef typename GeMatrix<MA>::ElementType  T;

    const Underscore<IndexType> _;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();
    const IndexType k = std::min(m, n);

    for (IndexType i=1; i<=k; ++i) {
//
//      Generate elementary reflector H(i) to annihilate A(i+1:m,i)
//
        larfg(m-i+1, A(i,i), A(_(std::min(i+1,m),m), i), tau(i));

        if (i<n) {
//
//          Apply H(i) to A(i:m,i+1:n) from the left
//
            const T Aii = A(i,i);
            A(i,i) = T(1);
            auto _work = work(_(1, n-i));
            larf(Left, A(_(i,m),i), tau(i), A(_(i,m), _(i+1,n)), _work);
            A(i,i) = Aii;
        }
    }
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename MA, typename VTAU, typename VWORK>
void
qr2_native(GeMatrix<MA> &A, DenseVector<VTAU> &tau, DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MA>::ElementType  T;

    const INTEGER    M = A.numRows();
    const INTEGER    N = A.numCols();
    const INTEGER    LDA = A.leadingDimension();
    INTEGER          INFO;

    if (IsSame<T, DOUBLE>::value) {
        LAPACK_IMPL(dgeqr2)(&M, &N, A.data(), &LDA,
                            tau.data(), work.data(),
                            &INFO);
        assert(INFO==0);
    } else {
        ASSERT(0);
    }
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename MA, typename VTAU, typename VWORK>
void
qr2(GeMatrix<MA> &A, DenseVector<VTAU> &tau, DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(tau.firstIndex()==1);
    ASSERT(work.firstIndex()==1);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();
    const IndexType k = std::min(m, n);

    ASSERT(tau.length()>=k);
    ASSERT(work.length()>=n);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MA>::NoView       _A      = A;
    typename DenseVector<VTAU>::NoView  _tau    = tau;
    typename DenseVector<VTAU>::NoView  _work   = work;
#   endif

//
//  Call implementation
//
    qr2_generic(A, tau, work);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    qr2_native(_A, _tau, _work);

    bool failed = false;
    if (! isIdentical(A, _A, " A", "_A")) {
        std::cerr << "CXXLAPACK:  A = " << A << std::endl;
        std::cerr << "F77LAPACK: _A = " << _A << std::endl;
        failed = true;
    }

    if (! isIdentical(tau, _tau, " tau", "_tau")) {
        std::cerr << "CXXLAPACK:  tau = " << tau << std::endl;
        std::cerr << "F77LAPACK: _tau = " << _tau << std::endl;
        failed = true;
    }

    if (! isIdentical(work, _work, " work", "_work")) {
        std::cerr << "CXXLAPACK:  work = " << work << std::endl;
        std::cerr << "F77LAPACK: _work = " << _work << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VTAU, typename VWORK>
void
qr2(MA &&A, VTAU &&tau, VWORK &&work)
{
    return qr2(A, tau, work);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_QR_QR2_TCC
