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
       SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_TF2_TCC
#define FLENS_LAPACK_IMPL_TF2_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename MA, typename VP>
typename GeMatrix<MA>::IndexType
tf2_impl(GeMatrix<MA> &A, DenseVector<VP> &piv)
{
    using lapack::lamch;
    using std::abs;
    using std::min;

    typedef typename GeMatrix<MA>::IndexType         IndexType;
    typedef typename GeMatrix<MA>::ElementType       T;
    typedef typename ComplexTrait<T>::PrimitiveType  PT;

    const Underscore<IndexType> _;

    const T  Zero(0), One(1);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    IndexType info = 0;

//
//  Quick return if possible
//
    if (m==0 || n==0) {
        return info;
    }
//
//     Compute machine safe minimum 
//
    const PT safeMin = lamch<PT>(SafeMin);

    for (IndexType j=1; j<=min(m,n); ++j) {
//
//      Row range of current submatrix A(j:M, j:N)
//
        const auto rows = _(j, m);
//
//      Row and column range of trailing submatrix A(j+1:M, j+1:N)
//
        const auto _rows = _(j+1, m);
        const auto _cols = _(j+1, n);
//
//      Find pivot and test for singularity.
//
        IndexType jp = j - 1 + blas::iamax(A(rows,j));
        piv(j) = jp;
        if (A(jp,j)!=Zero) {
//
//          Apply the interchange to columns 1:N.
//
            if (j!=jp) {
                blas::swap(A(j,_), A(jp,_));
            }
//
//          Compute elements J+1:M of J-th column.
//
            if (j<m) {
                if (abs(A(j,j))>=safeMin) {
                    blas::scal(One/A(j, j), A(_rows,j));
                } else {
                    for (IndexType i=1; i<=m-j; ++i) {
                        A(j+i,j) /= A(j,j);
                    }
                }
            }
        } else {
            if (info==0) {
                info = j;
            }
        }
//
//      Update trailing submatrix A(j+1:M, j+1:N)
//
        if (j<min(m,n)) {
            blas::ru(-One, A(_rows,j), A(j,_cols), A(_rows,_cols));
        }
    }

    return info;
}

} // namespace generic

//== interface for external lapack =============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename VP>
typename GeMatrix<MA>::IndexType
tf2_impl(GeMatrix<MA> &A, DenseVector<VP> &piv)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    return cxxlapack::getf2<IndexType>(A.numRows(), A.numCols(),
                                       A.data(), A.leadingDimension(),
                                       piv.data());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename MA, typename VP>
typename GeMatrix<MA>::IndexType
tf2(GeMatrix<MA> &A, DenseVector<VP> &piv)
{
    typedef typename GeMatrix<MA>::IndexType    IndexType;

//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT((piv.inc()>0 && piv.firstIndex()==1)
        || (piv.inc()<0 && piv.firstIndex()==A.numRows()));

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MA>::NoView       A_org      = A;
    typename DenseVector<VP>::NoView    piv_org    = piv;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::tf2_impl(A, piv);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results
//
    typename GeMatrix<MA>::NoView       A_generic      = A;
    typename DenseVector<VP>::NoView    piv_generic    = piv;
//
//  Restore output arguments
//
    A   = A_org;
    piv = piv_org;
//
//  Compare results
//
    IndexType _info = external::tf2_impl(A, piv);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(piv_generic, piv, "piv_generic", "piv")) {
        std::cerr << "CXXLAPACK: piv_generic = " << piv_generic << std::endl;
        std::cerr << "F77LAPACK: piv = " << piv << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, " info", "_info")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "A_org = " << A_org << std::endl;
        ASSERT(0);
    }
#   endif

    return info;
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VP>
typename MA::IndexType
tf2(MA &&A, VP &&piv)
{
    typedef typename MA::IndexType  IndexType;

    CHECKPOINT_ENTER;
    IndexType info =  tf2(A, piv);
    CHECKPOINT_LEAVE;

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_TF2_TCC
