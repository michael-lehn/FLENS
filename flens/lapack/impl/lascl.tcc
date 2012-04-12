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
       SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
 *
 *  -- LAPACK auxiliary routine (version 3.3.0) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2010
 */

#ifndef FLENS_LAPACK_IMPL_LASCL_TCC
#define FLENS_LAPACK_IMPL_LASCL_TCC 1

#include <cmath>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

////////
////////
////////
////////
////////
////////
////////
////////
////////        TODO: Completely redesign this crap!
////////
////////
////////
////////        ... but it works.
////////
////////
////////
////////
////////

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename IndexType, typename T, typename MA>
typename RestrictTo<IsSame<typename MA::ElementType, T>::value, void>::Type
lascl_impl(LASCL::Type   type,
           IndexType     kl,
           IndexType     ku,
           const T       &cFrom,
           const T       &cTo,
           MA            &A)
{
    using namespace LASCL;

    using std::abs;
    using std::isnan;
    using std::max;
    using std::min;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();
    const T Zero = 0;
    const T One = 1;

    if ((cFrom==Zero) || isnan(cFrom)) {
        ASSERT(0);
    } else if (isnan(cTo)) {
        ASSERT(0);
    } else if (type==SymmetricLowerBand) {
        ASSERT(A.numRows()==A.numCols());
    } else if (type==SymmetricUpperBand) {
        ASSERT(A.numRows()==A.numCols());
    }
//
//  Quick return if possible
//
    if ((n==0) || (m==0)) {
        return;
    }
//
//  Get machine parameters
//
    const T smallNum = lamch<T>(SafeMin);
    const T bigNum = T(1) / smallNum;
//
//  Make copies of cFrom, cTo
//
    T cFromC = cFrom;
    T cToC = cTo;

    bool done = false;

    do {
        T cFrom1 = cFromC*smallNum;
        T cTo1;
        T mul;
        if (cFrom1==cFromC) {
//          cFromC is an inf.  Multiply by a correctly signed zero for
//          finite cToC, or a NaN if cToC is infinite.
            mul = cToC / cFromC;
            done = true;
            cTo1 = cToC;
        } else {
            cTo1 = cToC / bigNum;
            if (cTo1==cToC) {
//              cToC is either 0 or an inf.  In both cases, cToC itself
//              serves as the correct multiplication factor.
                mul = cToC;
                done = true;
                cFromC = One;
            } else if (abs(cFrom1)>abs(cToC) && cToC!=Zero) {
                mul = smallNum;
                done = false;
                cFromC = cFrom1;
            } else if (abs(cTo1)>abs(cFromC)) {
                mul = bigNum;
                done = false;
                cToC = cTo1;
            } else {
                mul = cToC / cFromC;
                done = true;
            }
        }

        if (type==FullMatrix) {
//
//          Full matrix
//
            for (IndexType j=1; j<=n; ++j) {
                for (IndexType i=1; i<=m; ++i) {
                    A(i,j) *= mul;
                }
            }
        } else if (type==LowerTriangular) {
//
//          Lower triangular matrix
//
            for (IndexType j=1; j<=n; ++j) {
                for (IndexType i=j; i<=m; ++i) {
                    A(i,j) *= mul;
                }
            }
        } else if (type==UpperTriangular) {
//
//          Upper triangular matrix
//
            for (IndexType j=1; j<=n; ++j) {
                for (IndexType i=1; i<=min(j,m); ++i) {
                    A(i,j) *= mul;
                }
            }
        } else if (type==UpperHessenberg) {
//
//          Upper Hessenberg matrix
//
            for (IndexType j=1; j<=n; ++j) {
                for (IndexType i=1; i<=min(j+1,m); ++i) {
                    A(i,j) *= mul;
                }
            }
        } else if (type==SymmetricLowerBand) {
//
//          Lower half of a symmetric band matrix
//
            // TODO: this only works for the internal fullstorage of a
            //       band matrix not for the external element access
            //       of SbMatrix
            const IndexType k3 = kl + 1;
            const IndexType k4 = n + 1;
            for (IndexType j=1; j<=n; ++j) {
                for (IndexType i=1; i<=min(k3, k4-j); ++i) {
                    A(i,j) *= mul;
                }
            }
        } else if (type==SymmetricUpperBand) {
//
//          Upper half of a symmetric band matrix
//
            // TODO: this only works for the internal fullstorage of a
            //       band matrix not for the external element access
            //       of SbMatrix
            const IndexType k1 = ku + 2;
            const IndexType k3 = ku + 1;
            for (IndexType j=1; j<=n; ++j) {
                for (IndexType i=max(k1-j,IndexType(1)); i<=k3; ++i) {
                    A(i,j) *= mul;
                }
            }
        } else if (type==GeneralBand) {
//
//          Band matrix
//
            // TODO: this only works for the internal fullstorage of a
            //       band matrix not for the external element access
            //       of GeMatrix
            const IndexType k1 = kl + ku + 2;
            const IndexType k2 = kl + 1;
            const IndexType k3 = 2*kl + ku + 1;
            const IndexType k4 = kl + ku + 1 + m;
            for (IndexType j=1; j<=n; ++j) {
                for (IndexType i=max(k1-j,k2); i<=min(k3,k4-j); ++i) {
                    A(i,j) *= mul;
                }
            }
        } else {
            ASSERT(0);
        }
    } while (!done);
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename IndexType, typename T, typename MA>
void
lascl_impl(LASCL::Type   type,
           IndexType     kl,
           IndexType     ku,
           const T       &cFrom,
           const T       &cTo,
           MA            &A)
{
    cxxlapack::lascl<IndexType>(getF77Char(type),
                                kl,
                                ku,
                                cFrom,
                                cTo,
                                A.numRows(),
                                A.numCols(),
                                A.data(),
                                A.leadingDimension());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename IndexType, typename T, typename MA>
typename RestrictTo<IsSame<typename MA::ElementType, T>::value, void>::Type
lascl(LASCL::Type type, IndexType kl, IndexType ku,
      const T &cFrom, const T &cTo, MA &A)
{
    LAPACK_DEBUG_OUT("lascl");

#   ifdef CHECK_CXXLAPACK
    typename MA::NoView _A = A;
#   endif

    LAPACK_SELECT::lascl_impl(type, kl, ku, cFrom, cTo, A);

#   ifdef CHECK_CXXLAPACK
    external::lascl_impl(type, kl, ku, cFrom, cTo, _A);
    if (! isIdentical(A, _A, " A", "A_")) {
        std::cerr << "CXXLAPACK:  A = " << A << std::endl;
        std::cerr << "F77LAPACK: _A = " << _A << std::endl;
        ASSERT(0);
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename T, typename MA>
typename RestrictTo<IsSame<typename MA::ElementType, T>::value, void>::Type
lascl(LASCL::Type type, IndexType kl, IndexType ku,
      const T &cFrom, const T &cTo, MA &&A)
{
    lascl(type, kl, ku, cFrom, cTo, A);
}

//-- convert vector to matrix --------------------------------------------------
template <typename IndexType, typename T, typename VX>
void
lascl(LASCL::Type type, IndexType kl, IndexType ku,
      const T &cFrom, const T &cTo, DenseVector<VX> &x)
{
    using namespace LASCL;

    typedef typename DenseVector<VX>::ElementType   TX;
    typedef FullStorage<TX, ColMajor>               FS;

    typename GeMatrix<FS>::View  A(x.length(), 1, x);
    ASSERT(type==FullMatrix);
    lascl(type, kl, ku, cFrom, cTo, A);
}

//-- convert scalar to matrix --------------------------------------------------
template <typename IndexType, typename T, typename MA>
typename RestrictTo<IsSame<MA, T>::value, void>::Type
lascl(LASCL::Type   type,
      IndexType     kl,
      IndexType     ku,
      const T       &cFrom,
      const T       &cTo,
      MA            &scalar)
{
    using namespace LASCL;

    GeMatrixView<MA> A = typename GeMatrixView<MA>::Engine(1, 1, &scalar, 1);

    ASSERT(type==FullMatrix);
    lascl(type, kl, ku, cFrom, cTo, A);
}

template <typename IndexType, typename T, typename VX>
void
lascl(LASCL::Type type, IndexType kl, IndexType ku,
      const T &cFrom, const T &cTo, DenseVector<VX> &&x)
{
    lascl(type, kl, ku, cFrom, cTo, x);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_LASCL_TCC
