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

#ifndef FLENS_LAPACK_QR_ORM2R_TCC
#define FLENS_LAPACK_QR_ORM2R_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VTAU, typename MC, typename VWORK>
void
orm2r(Side side, Transpose trans, MA &&A, const VTAU &tau, MC &&C,
      VWORK &&work)
{
    orm2r(side, trans, A, tau, C, work);
}

//-- ormqr ---------------------------------------------------------------------
template <typename MA, typename VTAU, typename MC, typename VWORK>
void
orm2r(Side side, Transpose trans, GeMatrix<MA> &A,
      const DenseVector<VTAU> &tau, GeMatrix<MC> &C,
      DenseVector<VWORK> &work)
{
#   ifndef NDEBUG
    if ((side==Left) && (work.length()<C.numCols())) {
        ASSERT(0);
    }
    if ((side==Right) && (work.length()<C.numRows())) {
        ASSERT(0);
    }
#   endif

    typedef typename GeMatrix<MC>::IndexType    IndexType;
    typedef typename GeMatrix<MC>::ElementType  T;

    typedef Range<IndexType>    Range;
    const Underscore<IndexType> _;

    IndexType m = C.numRows();
    IndexType n = C.numCols();
    IndexType k = A.numCols();

//
//  nq is the order of Q
//
    const IndexType nq = (side==Left) ? m : n;
    
    ASSERT(A.numRows()>=nq);
    ASSERT(k<=nq);

//
//  Quick return if possible
//
    if ((m==0) || (n==0) || (k==0)) {
        return;
    }

    IndexType iBeg, iEnd, iInc;
    if (((side==Left) && (trans!=NoTrans))
     || ((side==Right) && (trans==NoTrans)))
    {
        iBeg = 1;
        iEnd = k+1;
        iInc = 1;
    } else {
        iBeg = k;
        iEnd = 0;
        iInc = -1;
    }

    Range rows = _(1,m);
    Range cols = _(1,n);

    for (IndexType i=iBeg; i!=iEnd; i+=iInc) {
        if (side==Left) {
//
//          H(i) is applied to C(i:m,1:n)
//
            rows = _(i,m);
        } else {
//
//          H(i) is applied to C(1:m,i:n)
//
            cols = _(i,n);
        }
//
//      Apply H(i)
//
        const T Aii = A(i,i);
        A(i,i) = T(1);
        larf(side, A(_(i,m), i), tau(i), C(rows, cols), work);
        A(i,i) = Aii;
    }
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_QR_ORM2R_TCC
