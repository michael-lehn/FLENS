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

/* Baesed on
 *
      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )

      SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
      SUBROUTINE ZTRTRI( UPLO, DIAG, N, A, LDA, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_GE_TRI_TCC
#define FLENS_LAPACK_GE_TRI_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- (ge)tri [real and complex variant] ----------------------------------------

template <typename MA, typename VP, typename VWORK>
typename GeMatrix<MA>::IndexType
tri_impl(GeMatrix<MA>            &A,
         const DenseVector<VP>   &piv,
         DenseVector<VWORK>      &work)
{
    using std::max;

    typedef typename GeMatrix<MA>::ElementType ElementType;
    typedef typename GeMatrix<MA>::IndexType   IndexType;

    const ElementType Zero(0), One(1);
    const IndexType n = A.numRows();
    const Underscore<IndexType> _;

    IndexType info = 0;
    IndexType nb = ilaenv<ElementType>(1, "GETRI", "", n);

    const IndexType lWorkOpt = n*nb;

    if (work.length()==0) {
        work.resize(max(lWorkOpt, IndexType(1)));
    }
    work(1) = lWorkOpt;

//
//  Quick return if possible
//
    if (n==0) {
        return info;
    }
//
//  Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,
//  and the inverse is not computed.
//
    info = tri(A.upper());
    if (info>0) {
        return info;
    }

    IndexType nbMin = 2;
    const IndexType lWork  = work.length();;
    const IndexType ldWork = n;

    IndexType iws;

    if (nb>1 && nb<n) {
        iws = max(ldWork*nb, IndexType(1));
        if (lWork<iws) {
            nb = lWork / ldWork;
            nbMin = max(2, ilaenv<ElementType>(2, "GETRI", "", n));
        }
    } else {
        iws = n;
    }

    GeMatrixView<ElementType> Work(n, nb, work);
//
//  Solve the equation inv(A)*L = inv(U) for inv(A).
//
    if (nb<nbMin || nb>=n) {
//
//      Use unblocked code.
//
        for (IndexType j=n; j>=1; --j) {
//
//          Copy current column of L to WORK and replace with zeros.
//
            work(_(j+1,n)) = A(_(j+1,n),j);
            A(_(j+1,n),j)  = Zero;;
//
//          Compute current column of inv(A).
//
            if (j<n) {
                blas::mv(NoTrans, -One,
                         A(_,_(j+1,n)), work(_(j+1,n)),
                         One,
                         A(_,j));
            }
        }
    } else {
//
//      Use blocked code.
//
        const IndexType nn = ((n-1)/nb)*nb + 1;
        for (IndexType j=nn; j>=1; j-=nb) {
            const IndexType jb = min(nb, n-j+1);
//
//          Copy current block column of L to WORK and replace with
//          zeros.
//
            for (IndexType jj=j, JJ=1; jj<=j+jb-1; ++jj, ++JJ) {
                Work(_(jj+1,n),JJ) = A(_(jj+1,n),jj);
                A(_(jj+1,n),jj)    = Zero;
            }
//
//          Compute current block column of inv(A).
//
            if (j+jb<=n) {
                blas::mm(NoTrans, NoTrans,
                         -One,
                         A(_,_(j+jb,n)),
                         Work(_(j+jb,n),_(1,jb)),
                         One,
                         A(_,_(j,j+jb-1)));
            }
            blas::sm(Right, NoTrans,
                     One, Work(_(j,j+jb-1),_(1,jb)).lowerUnit(),
                     A(_,_(j,j+jb-1)));
        }
    }
//
//  Apply column interchanges.
//
    for (IndexType j=n-1; j>=1; --j) {
        const IndexType jp = piv(j);
        if (jp!=j) {
            blas::swap(A(_,j), A(_,jp));
        }
    }

    work(1) = iws;
    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (ge)tri [real and complex variant] ----------------------------------------

template <typename MA, typename VP, typename VWORK>
typename GeMatrix<MA>::IndexType
tri_impl(GeMatrix<MA>             &A,
         const DenseVector<VP>    &piv,
         DenseVector<VWORK>       &work)
{
    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    IndexType info;

    if (work.length()==0) {
        ElementType  WORK;
        IndexType    LWORK = -1;

        info = cxxlapack::getri<IndexType>(A.numRows(),
                                           A.data(),
                                           A.leadingDimension(),
                                           piv.data(),
                                           &WORK,
                                           LWORK);
        work.resize(IndexType(cxxblas::real(WORK)));
    }

    info = cxxlapack::getri<IndexType>(A.numRows(),
                                       A.data(),
                                       A.leadingDimension(),
                                       piv.data(),
                                       work.data(),
                                       work.length());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//-- (ge)tri [real and complex variant] ----------------------------------------

template <typename MA, typename VPIV, typename VWORK>
typename RestrictTo<(IsRealGeMatrix<MA>::value
                  && IsIntegerDenseVector<VPIV>::value
                  && IsRealDenseVector<VWORK>::value)
               ||   (IsComplexGeMatrix<MA>::value
                  && IsIntegerDenseVector<VPIV>::value
                  && IsComplexDenseVector<VWORK>::value),
         typename RemoveRef<MA>::Type::IndexType>::Type
tri(MA          &&A,
    const VPIV  &piv,
    VWORK       &&work)
{
    using std::max;

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<VWORK>::Type VectorWork;
#   endif

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());

    const IndexType n = A.numRows();

    ASSERT(piv.firstIndex()==1);
    ASSERT(piv.length()==n);

    const bool lQuery = (work.length()==0);
    ASSERT(lQuery || work.length()>=max(IndexType(1),n));
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView        A_org    = A;
    typename VectorWork::NoView     work_org = work;
#   endif

//
//  Call implementation
//
    const IndexType info = LAPACK_SELECT::tri_impl(A, piv, work);

//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView        A_generic    = A;
    typename VectorWork::NoView     work_generic = work;

    A    = A_org;
    work = work_org;

    const IndexType _info = external::tri_impl(A, piv, work);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, " info", "_info")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    } else {
        // std::cerr << "passed: (ge)tri.tcc" << std::endl;
    }
#   endif

    return info;
}

//-- (ge)tri [real/complex variant with temporary workspace] -------------------
template <typename MA, typename VPIV>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
tri(MA          &&A,
    const VPIV  &piv)
{
    typedef typename RemoveRef<MA>::Type::Vector WorkVector;

    WorkVector  work;
    return tri(A, piv, work);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GE_TRI_TCC
