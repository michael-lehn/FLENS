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

#ifndef FLENS_LAPACK_TR_TRI_TCC
#define FLENS_LAPACK_TR_TRI_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- (tr)tri [real and complex variant] ----------------------------------------

template <typename MA>
typename GeMatrix<MA>::IndexType
tri_impl(TrMatrix<MA> &A)
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

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (tr)tri [real and complex variant] ----------------------------------------

template <typename MA>
typename GeMatrix<MA>::IndexType
tri_impl(TrMatrix<MA> &A)
{
    typedef typename TrMatrix<MA>::IndexType  IndexType;

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

//-- (tr)tri [real and complex variant] ----------------------------------------

template <typename MA>
typename RestrictTo<IsRealTrMatrix<MA>::value
                 || IsComplexTrMatrix<MA>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
tri(MA &&A)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

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
    typename MatrixA::NoView   A_org = A;
#   endif

//
//  Call implementation
//
    const IndexType info = LAPACK_SELECT::tri_impl(A);

//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView   A_generic = A;

    A = A_org;

    const IndexType _info = external::tri_impl(A);

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

} } // namespace lapack, flens

#endif // FLENS_LAPACK_TR_TRI_TCC
