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

/* Baesed on
 *
      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )

      SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_GESV_TRI_TCC
#define FLENS_LAPACK_GESV_TRI_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
//-- (ge)tri
template <typename MA, typename VP, typename VWORK>
typename GeMatrix<MA>::IndexType
tri_generic(GeMatrix<MA>            &A,
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

//-- (tr)tri
template <typename MA>
typename GeMatrix<MA>::IndexType
tri_generic(TrMatrix<MA> &A)
{
    using std::min;
    using cxxblas::getF77BlasChar;

    typedef typename TrMatrix<MA>::ElementType ElementType;
    typedef typename TrMatrix<MA>::IndexType   IndexType;

    const Underscore<IndexType> _;

    const IndexType n = A.dim();
    const bool upper  = (A.upLo()==Upper);
    const bool noUnit = (A.diag()==NonUnit);

    const ElementType  Zero(0), One(1);

    IndexType info = 0;
//
//  Quick return if possible
//
    if (n==0) {
        return info;
    }
//
//  Check for singularity if non-unit.
//
    if (noUnit) {
        for (info=1; info<=n; ++info) {
            if (A(info,info)==Zero) {
                return info;
            }
        }
        info = 0;
    }

//
//  Determine the block size for this environment.  //
    const char upLoDiag[2] = { getF77BlasChar(A.upLo()),
                               getF77BlasChar(A.diag()) };
    const IndexType nb = ilaenv<ElementType>(1, "TRTRI", upLoDiag, n);

    if (nb<=1 || nb>=n) {
//
//      Use unblocked code
//
        info = ti2(A);
    } else {
//
//      Use blocked code
//

        if (upper) {
//
//          Compute inverse of upper triangular matrix
//
            for (IndexType j=1; j<=n; j+=nb) {
                const IndexType jb = min(nb, n-j+1);
//
//              Compute rows 1:j-1 of current block column
//
                const auto range1 = _(1,j-1);
                const auto range2 = _(j,j+jb-1);

                const auto U11 = (noUnit) ? A(range1, range1).upper()
                                          : A(range1, range1).upperUnit();
                auto U22 = (noUnit) ? A(range2, range2).upper()
                                    : A(range2, range2).upperUnit();
                auto U12 = A(range1, range2);

                blas::mm(Left, NoTrans, One, U11, U12);
                blas::sm(Right, NoTrans, -One, U22, U12);
//
//              Compute inverse of current diagonal block
//
                info = ti2(U22);
            }
        } else {
//
//          Compute inverse of lower triangular matrix
//
            const IndexType nn = ((n-1)/nb)*nb + 1;
            for (IndexType j=nn; j>=1; j-=nb) {
                const IndexType jb = min(nb, n-j+1);

                const auto range1 = _(j,j+jb-1);
                auto L11 = (noUnit) ? A(range1,range1).lower()
                                    : A(range1,range1).lowerUnit();

                if (j+jb<=n) {
//
//                  Compute rows j+jb:n of current block column
//
                    const auto range2 = _(j+jb,n);

                    const auto L22 = (noUnit) ? A(range2,range2).lower()
                                              : A(range2,range2).lowerUnit();

                    auto L21 = A(range2, range1);

                    blas::mm(Left, NoTrans, One, L22, L21);
                    blas::sm(Right, NoTrans, -One, L11, L21);
                }
//
//              Compute inverse of current diagonal block
//
                info = ti2(L11);
            }
        }
    }
    return info;
}

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (ge)tri
template <typename MA, typename VP, typename VWORK>
typename GeMatrix<MA>::IndexType
tri(GeMatrix<MA>             &A,
    const DenseVector<VP>    &piv,
    DenseVector<VWORK>       &work)
{
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    IndexType info = cxxlapack::getri<IndexType>(A.numRows(),
                                                 A.data(),
                                                 A.leadingDimension(),
                                                 piv.data(),
                                                 work.data(),
                                                 work.length());
    ASSERT(info>=0);
    return info;
}

//-- (tr)tri
template <typename MA>
typename GeMatrix<MA>::IndexType
tri(TrMatrix<MA> &A)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::trtri<IndexType>(getF77Char(A.upLo()),
                                                 getF77Char(A.diag()),
                                                 A.dim(),
                                                 A.data(),
                                                 A.leadingDimension());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//-- (ge)tri
template <typename MA, typename VP, typename VWORK>
typename GeMatrix<MA>::IndexType
tri(GeMatrix<MA> &A, const DenseVector<VP> &piv, DenseVector<VWORK> &work)
{
    using std::max;

    typedef typename GeMatrix<MA>::IndexType IndexType;
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
    ASSERT(lQuery || work.length()>=n);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename GeMatrix<MA>::NoView        A_org    = A;
    typename DenseVector<VWORK>::NoView  work_org = work;
#   endif

//
//  Call implementation
//
    const IndexType info = tri_generic(A, piv, work);

//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    typename GeMatrix<MA>::NoView        A_generic    = A;
    typename DenseVector<VWORK>::NoView  work_generic = work;

    A    = A_org;
    work = work_org;

    const IndexType _info = cxxlapack::tri(A, piv, work);

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

//-- (tr)tri
template <typename MA>
typename GeMatrix<MA>::IndexType
tri(TrMatrix<MA> &A)
{
    typedef typename GeMatrix<MA>::IndexType IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename TrMatrix<MA>::NoView   A_org = A;
#   endif

//
//  Call implementation
//
    const IndexType info = tri_generic(A);

//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    typename TrMatrix<MA>::NoView   A_generic = A;

    A = A_org;

    const IndexType _info = cxxlapack::tri(A);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
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
        // std::cerr << "passed: (tr)tri.tcc" << std::endl;
    }
#   endif

    return info;
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VP, typename VWORK>
typename MA::IndexType
tri(MA &&A, const VP &&piv, VWORK &&work)
{
    typedef typename MA::IndexType  IndexType;

    CHECKPOINT_ENTER;
    IndexType info = tri(A, piv, work);
    CHECKPOINT_LEAVE;

    return info;
}

template <typename MA>
typename MA::IndexType
tri(MA &&A)
{
    typedef typename MA::IndexType  IndexType;

    CHECKPOINT_ENTER;
    IndexType info = tri(A);
    CHECKPOINT_LEAVE;

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GESV_TRI_TCC
