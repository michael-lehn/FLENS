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
       SUBROUTINE DGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,
     $                   INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_EIG_EQU_TCC
#define FLENS_LAPACK_EIG_EQU_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
template <typename MA, typename VR, typename VC,
          typename ROWCOND, typename COLCOND,
          typename AMAX>
typename GeMatrix<MA>::IndexType
equ_generic(const GeMatrix<MA>  &A,
            DenseVector<VR>     &r,
            DenseVector<VC>     &c,
            ROWCOND             &rowCond,
            COLCOND             &colCond,
            AMAX                &amax)
{
    using std::abs;
    using std::max;
    using std::min;

    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const T  Zero(0), One(1);

    const Underscore<IndexType>  _;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    IndexType info = 0;
//
//  Quick return if possible
//
    if (m==0 || n==0) {
        rowCond = One;
        colCond = One;
        amax = Zero;
        return info;
    }
//
//  Get machine constants.
//
    const T smallNum = lamch<T>(SafeMin);
    const T bigNum = One / smallNum;
//
//  Compute row scale factors.
//
    r = Zero;
//
//  Find the maximum element in each row.
//
    for (IndexType j=1; j<=n; ++j) {
        for (IndexType i=1; i<=m; ++i) {
            r(i) = max(r(i),abs(A(i,j)));
        }
    }
//
//  Find the maximum and minimum scale factors.
//
    T rcMin = bigNum;
    T rcMax = Zero;
    for (IndexType i=1; i<=m; ++i) {
        rcMax = max(rcMax, r(i));
        rcMin = min(rcMin, r(i));
    }
    amax = rcMax;

    if (rcMin==Zero) {
//
//      Find the first zero scale factor and return an error code.
//
        for (IndexType i=1; i<=m; ++i) {
            if (r(i)==Zero) {
                info = i;
                return info;
            }
        }
    } else {
//
//      Invert the scale factors.
//
        for (IndexType i=1; i<=m; ++i) {
            r(i) = One / min(max(r(i),smallNum), bigNum);
        }
//
//      Compute ROWCND = min(R(I)) / max(R(I))
//
        rowCond = max(rcMin,smallNum) / min(rcMax,bigNum);
    }

//  Compute column scale factors
//
    c = Zero;
//
//  Find the maximum element in each column,
//  assuming the row scaling computed above.
//
    for (IndexType j=1; j<=n; ++j) {
        for (IndexType i=1; i<=m; ++i) {
            c(j) = max(c(j),abs(A(i,j))*r(i));
        }
    }
//
//  Find the maximum and minimum scale factors.
//
    rcMin = bigNum;
    rcMax = Zero;
    for (IndexType j=1; j<=n; ++j) {
        rcMin = min(rcMin, c(j));
        rcMax = max(rcMax, c(j));
    }

    if (rcMin==Zero) {
//
//      Find the first zero scale factor and return an error code.
//
        for (IndexType j=1; j<=n; ++j) {
            if (c(j)==Zero) {
                info = m+j;
                return info;
            }
        }
    } else {
//
//      Invert the scale factors.
//
        for (IndexType j=1; j<=n; ++j) {
            c(j) = One / min(max(c(j),smallNum), bigNum);
        }
//
//      Compute COLCND = min(C(J)) / max(C(J))
//
        colCond = max(rcMin,smallNum) / min(rcMax,bigNum);
    }
    return info;
}


//== interface for native lapack ===============================================


#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename VR, typename VC,
          typename ROWCOND, typename COLCOND,
          typename AMAX>
typename GeMatrix<MA>::IndexType
equ(const GeMatrix<MA>  &A,
    DenseVector<VR>     &r,
    DenseVector<VC>     &c,
    ROWCOND             &rowCond,
    COLCOND             &colCond,
    AMAX                &amax)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::geequ<IndexType>(A.numRows(),
                                                 A.numCols(),
                                                 A.data(),
                                                 A.leadingDimension(),
                                                 r.data(),
                                                 c.data(),
                                                 &rowCond,
                                                 &colCond,
                                                 &amax,
                                                 &info);
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename MA, typename VR, typename VC,
          typename ROWCOND, typename COLCOND,
          typename AMAX>
typename GeMatrix<MA>::IndexType
equ(const GeMatrix<MA>  &A,
    DenseVector<VR>     &r,
    DenseVector<VC>     &c,
    ROWCOND             &rowCond,
    COLCOND             &colCond,
    AMAX                &amax)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    ASSERT(r.firstIndex()==1);
    ASSERT(r.length()==m);

    ASSERT(c.firstIndex()==1);
    ASSERT(c.length()==n);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename DenseVector<VR>::NoView  r_org       = r;
    typename DenseVector<VC>::NoView  c_org       = c;
    ROWCOND                           rowCond_org = rowCond;
    COLCOND                           colCond_org = colCond;
    AMAX                              amax_org    = amax;
#   endif

//
//  Call implementation
//
    const IndexType info = equ_generic(A, r, c, rowCond, colCond, amax);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename DenseVector<VR>::NoView  r_generic       = r;
    typename DenseVector<VC>::NoView  c_generic       = c;
    ROWCOND                           rowCond_generic = rowCond;
    COLCOND                           colCond_generic = colCond;
    AMAX                              amax_generic    = amax;
//
//  restore output arguments
//
    r = r_org;
    c = c_org;
    rowCond = rowCond_org;
    colCond = colCond_org;
    amax = amax_org;

//
//  Compare generic results with results from the native implementation
//
    const IndexType info_ = external::equ(A, r, c, rowCond, colCond, amax);

    bool failed = false;
    if (! isIdentical(r_generic, r, "r_generic", "r")) {
        std::cerr << "CXXLAPACK: r_generic = " << r_generic << std::endl;
        std::cerr << "F77LAPACK: r = " << r << std::endl;
        failed = true;
    }

    if (! isIdentical(c_generic, c, "c_generic", "c")) {
        std::cerr << "CXXLAPACK: c_generic = " << c_generic << std::endl;
        std::cerr << "F77LAPACK: c = " << c << std::endl;
        failed = true;
    }

    if (! isIdentical(rowCond_generic, rowCond, "rowCond_generic", "rowCond")) {
        std::cerr << "CXXLAPACK: rowCond_generic = "
                  << rowCond_generic << std::endl;
        std::cerr << "F77LAPACK: rowCond = "
                  << rowCond << std::endl;
        failed = true;
    }

    if (! isIdentical(colCond_generic, colCond, "colCond_generic", "colCond")) {
        std::cerr << "CXXLAPACK: colCond_generic = "
                  << colCond_generic << std::endl;
        std::cerr << "F77LAPACK: colCond = "
                  << colCond << std::endl;
        failed = true;
    }

    if (! isIdentical(amax_generic, amax, "amax_generic", "amax")) {
        std::cerr << "CXXLAPACK: amax_generic = "
                  << amax_generic << std::endl;
        std::cerr << "F77LAPACK: amax = " << amax << std::endl;
        failed = true;
    }

    if (info!=info_) {
        std::cerr << "CXXLAPACK: info = " << info << std::endl;
        std::cerr << "F77LAPACK: info_ = " << info_ << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "A = " << A << std::endl;
        std::cerr << "r = " << r << std::endl;
        std::cerr << "c = " << c << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: equ.tcc" << std::endl;
    }
#   endif

    return info;
}


//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VR, typename VC,
          typename ROWCOND, typename COLCOND,
          typename AMAX>
typename MA::IndexType
equ(const MA  &A,
    VR        &&r,
    VC        &&c,
    ROWCOND   &&rowCond,
    COLCOND   &&colCond,
    AMAX      &&amax)
{
    CHECKPOINT_ENTER;
    const auto info = equ(A, r, c, rowCond, colCond, amax);
    CHECKPOINT_LEAVE;

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_EQU_TCC
