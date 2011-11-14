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
       DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_AUX_LANGE_TCC
#define FLENS_LAPACK_AUX_LANGE_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename MA, typename VWORK>
typename GeMatrix<MA>::ElementType
lange_generic(Norm norm, const GeMatrix<MA> &A, DenseVector<VWORK> &work)
{
    using std::abs;
    using std::max;
    using std::min;
    using std::sqrt;

    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const Underscore<IndexType> _;
    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    if (min(m,n)==T(0)) {
        return 0;
    } else if (norm==MaximumNorm) {
//
//      Find max(abs(A(i,j))).
//
        T value = 0;

        for (IndexType j=1; j<=n; ++j) {
            for (IndexType i=1; i<=m; ++i) {
                value = max(value, abs(A(i,j)));
            }
        }
        return value;
    } else if (norm==OneNorm) {
//
//      Find norm1(A).
//
        T value = 0;

        for (IndexType j=1; j<=n; ++j) {
            T sum = 0;
            for (IndexType i=1; i<=m; ++i) {
                sum += abs(A(i,j));
            }
            value = max(value, sum);
        }
        return value;
    } else if (norm==InfinityNorm) {
//
//      Find normI(A).
//
        for (IndexType i=1; i<=m; ++i) {
            work(i) = 0;
        }
        for (IndexType j=1; j<=n; ++j) {
            for (IndexType i=1; i<=m; ++i) {
                work(i) += abs(A(i,j));
            }
        }
        T value = 0;
        for (IndexType i=1; i<=m; ++i) {
            value = max(value, work(i));
        }
        return value;
    } else if (norm==FrobeniusNorm) {
//
//      Find normF(A).
//
        T scale = 0;
        T sum = 1;
        for (IndexType j=1; j<=n; ++j) {
            lassq(A(_,j), scale, sum);
        }
        return scale*sqrt(sum);
    }
    ASSERT(0);
    return 0;
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename MA, typename VWORK>
typename GeMatrix<MA>::ElementType
lange_native(Norm norm, const GeMatrix<MA> &A, DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MA>::ElementType  T;

    const char      NORM    = getF77LapackChar<Norm>(norm);
    const INTEGER   M       = A.numRows();
    const INTEGER   N       = A.numCols();
    const INTEGER   LDA     = A.leadingDimension();

    if (IsSame<T, DOUBLE>::value) {
        return LAPACK_IMPL(dlange)(&NORM,
                                   &M,
                                   &N,
                                   A.data(),
                                   &LDA,
                                   work.data());
    } else {
        ASSERT(0);
    }
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename MA>
typename GeMatrix<MA>::ElementType
lange(Norm norm, const GeMatrix<MA> &A)
{
    ASSERT(norm!=InfinityNorm);

    typedef typename GeMatrix<MA>::ElementType  T;
    DenseVector<Array<T> >  dummy;
    return lange(norm, A, dummy);
}

template <typename MA, typename VWORK>
typename GeMatrix<MA>::ElementType
lange(Norm norm, const GeMatrix<MA> &A, DenseVector<VWORK> &work)
{
    LAPACK_DEBUG_OUT("lange");

    typedef typename GeMatrix<MA>::ElementType T;
//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(norm!=InfinityNorm || work.length()>=A.numRows());
#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename DenseVector<VWORK>::NoView _work = work;
#   endif

//
//  Call implementation
//
    T result = lange_generic(norm, A, work);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    if (_work.length()==0) {
        _work.resize(work.length());
    }

    T _result = lange_native(norm, A, _work);

    bool failed = false;
    if (! isIdentical(work, _work, " work", "_work")) {
        std::cerr << "CXXLAPACK:  work = " << work << std::endl;
        std::cerr << "F77LAPACK: _work = " << _work << std::endl;
        failed = true;
    }

    if (! isIdentical(result, _result, " result", "_result")) {
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif

    return result;
}

//-- forwarding ----------------------------------------------------------------
template <typename MA>
typename MA::ElementType
lange(Norm norm, const MA &A)
{
    typedef typename MA::ElementType  T;

    CHECKPOINT_ENTER;
    const T result = lange(norm, A);
    CHECKPOINT_LEAVE;

    return result;
}

template <typename MA, typename VWORK>
typename MA::ElementType
lange(Norm norm, const MA &A, VWORK &&work)
{
    typedef typename MA::ElementType  T;

    CHECKPOINT_ENTER;
    const T result = lange(norm, A, work);
    CHECKPOINT_LEAVE;

    return result;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_LANGE_TCC
