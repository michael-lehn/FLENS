/*
 *   Copyright (c) 2012, Michael Lehn, Klaus Pototzky
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
       SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
       SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_GE_TRF_TCC
#define FLENS_LAPACK_GE_TRF_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- (ge)trf [real and complex variant] ----------------------------------------

template <typename MA, typename VP>
typename GeMatrix<MA>::IndexType
trf_impl(GeMatrix<MA> &A, DenseVector<VP> &piv)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const Underscore<IndexType> _;

    const IndexType m  = A.numRows();
    const IndexType n  = A.numCols();
    const IndexType mn = min(m,n);

    IndexType info = 0;
//
//  Quick return if possible
//
    if ((m==0) || (n==0)) {
        return info;
    }
//
//  Determine the block size for this environment.
//
    IndexType bs = ilaenv<T>(1, "GETRF", "", m, n, -1, -1);

    if ((bs<=1) || (bs>=mn)) {
//
//      Use unblocked code.
//
        info = tf2(A, piv);
    } else {
//
//      Use blocked code.
//
        for (IndexType j=1; j<=mn; j+=bs) {
            IndexType jb = min(mn-j+1, bs);
//
//          Row and column partitioning of A
//
            const auto rows1    = _(    j, j+jb-1);
            const auto rows2    = _( j+jb,      m);

            const auto rows12   = _(    j,      m);

            const auto cols0    = _(    1,    j-1);
            const auto cols1    = _(    j, j+jb-1);
            const auto cols2    = _( j+jb,      n);
//
//          Factor diagonal and subdiagonal blocks and test for exact
//          singularity.
//
            IndexType _info = tf2(A(rows12, cols1), piv(_(j,mn)));
//
//          Adjust INFO and the pivot indices.
//
            if ((info==0) && (_info>0)) {
                info = _info + j - 1;
            }
            for (IndexType i=j; i<=min(m,j+jb-1); ++i) {
                piv(i) += j-1;
            }
//
//          Apply interchanges to columns 1:J-1.
//
            laswp(A(_,cols0), piv(rows1, j));

            if (j+jb<=n) {
//
//              Apply interchanges to columns J+JB:N.
//
                laswp(A(_,cols2), piv(rows1, j));
//
//              Compute block row of U.
//
                blas::sm(Left, NoTrans, T(1),
                         A(rows1, cols1).lowerUnit(),
                         A(rows1, cols2));

                if (j+jb<=m) {
//
//                  Update trailing submatrix.
//
                    blas::mm(NoTrans, NoTrans,
                             T(-1), A(rows2, cols1), A(rows1, cols2),
                             T( 1), A(rows2, cols2));
                }
            }
        }
    }
    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (ge)trf [real and complex variant] ----------------------------------------

template <typename MA, typename VP>
typename GeMatrix<MA>::IndexType
trf_impl(GeMatrix<MA> &A, DenseVector<VP> &piv)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    return cxxlapack::getrf<IndexType>(A.numRows(), A.numCols(),
                                       A.data(), A.leadingDimension(),
                                       piv.data());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//-- (ge)trf [real and complex variant] ----------------------------------------

template <typename MA, typename VPIV>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
trf(MA &&A, VPIV &&piv)
{
    using std::min;

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename RemoveRef<VPIV>::Type  VectorPiv;

//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT((piv.inc()>0 && piv.firstIndex()==1)
        || (piv.inc()<0 && piv.firstIndex()==A.numRows()));

    const IndexType mn = min(A.numRows(), A.numCols());

    if (piv.length()!=mn) {
        std::cerr << "A.numRows() = " << A.numRows() << std::endl;
        std::cerr << "A.numCols() = " << A.numCols() << std::endl;
        std::cerr << "mn = " << mn << std::endl;
        std::cerr << "piv.length() = " << piv.length() << std::endl;
    }
    if (piv.length()==0) {
        piv.resize(mn);
    }
    if (piv.length()!=mn) {
        std::cerr << "piv.length() = " << piv.length() << std::endl;
    }
    ASSERT(piv.length()==mn);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename MatrixA::NoView    _A      = A;
    typename VectorPiv::NoView  _piv    = piv;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::trf_impl(A, piv);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    IndexType _info = external::trf_impl(_A, _piv);

    bool failed = false;
    if (! isIdentical(A, _A, " A", "_A")) {
        std::cerr << "CXXLAPACK:  A = " << A << std::endl;
        std::cerr << "F77LAPACK: _A = " << _A << std::endl;
        failed = true;
    }

    if (! isIdentical(piv, _piv, " piv", "_piv")) {
        std::cerr << "CXXLAPACK:  piv = " << piv << std::endl;
        std::cerr << "F77LAPACK: _piv = " << _piv << std::endl;
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

#endif // FLENS_LAPACK_GE_TRF_TCC
