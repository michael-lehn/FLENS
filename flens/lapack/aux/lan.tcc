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

#ifndef FLENS_LAPACK_AUX_LAN_TCC
#define FLENS_LAPACK_AUX_LAN_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename MA, typename VWORK>
typename GeMatrix<MA>::ElementType
lan_generic(Norm norm, const GeMatrix<MA> &A, DenseVector<VWORK> &work)
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

    const T  Zero(0), One(1);

    if (min(m,n)==T(0)) {
        return Zero;
    } else if (norm==MaximumNorm) {
//
//      Find max(abs(A(i,j))).
//
        T value = Zero;

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
        T value = Zero;

        for (IndexType j=1; j<=n; ++j) {
            T sum = Zero;
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
            work(i) = Zero;
        }
        for (IndexType j=1; j<=n; ++j) {
            for (IndexType i=1; i<=m; ++i) {
                work(i) += abs(A(i,j));
            }
        }
        T value = Zero;
        for (IndexType i=1; i<=m; ++i) {
            value = max(value, work(i));
        }
        return value;
    } else if (norm==FrobeniusNorm) {
//
//      Find normF(A).
//
        T scale = Zero;
        T sum = One;
        for (IndexType j=1; j<=n; ++j) {
            lassq(A(_,j), scale, sum);
        }
        return scale*sqrt(sum);
    }
    ASSERT(0);
    return Zero;
}

template <typename MA, typename VWORK>
typename TrMatrix<MA>::ElementType
lan_generic(Norm norm, const TrMatrix<MA> &A, DenseVector<VWORK> &work)
{
    using std::abs;
    using std::max;
    using std::min;
    using std::sqrt;

    typedef typename TrMatrix<MA>::ElementType  T;
    typedef typename TrMatrix<MA>::IndexType    IndexType;

    const Underscore<IndexType> _;
    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    const T  Zero(0), One(1);

    T value = Zero;

    if (min(m,n)==0) {
        return value;
    } else if (norm==MaximumNorm) {
//
//      Find max(abs(A(i,j))).
//
        if (A.diag()==Unit) {
            value = One;
            if (A.upLo()==Upper) {
                for (IndexType j=1; j<=n; ++j) {
                    for (IndexType i=1; i<=min(m,j-1); ++i) {
                        value = max(value, abs(A(i,j)));
                    }
                }
            } else {
                for (IndexType j=1; j<=n; ++j) {
                    for (IndexType i=j+1; i<=m; ++i) {
                        value = max(value, abs(A(i,j)));
                    }
                }
            }
        } else {
            value = Zero;
            if (A.upLo()==Upper) {
                for (IndexType j=1; j<=n; ++j) {
                    for (IndexType i=1; i<=min(m,j); ++i) {
                        value = max(value, abs(A(i,j)));
                    }
                }
            } else {
                for (IndexType j=1; j<=n; ++j) {
                    for (IndexType i=j; i<=m; ++i) {
                        value = max(value, abs(A(i,j)));
                    }
                }
            }
        }
    } else if (norm==OneNorm) {
//
//      Find norm1(A).
//
        value = Zero;
        const bool unitDiag = (A.diag()==Unit);
        if (A.upLo()==Upper) {
            for (IndexType j=1; j<=n; ++j) {
                T sum;
                if (unitDiag && (j<=m)) {
                    sum = One;
                    for (IndexType i=1; i<=j-1; ++i) {
                        sum += abs(A(i,j));
                    }
                } else {
                    sum = Zero;
                    for (IndexType i=1; i<=min(m,j); ++i) {
                        sum += abs(A(i,j));
                    }
                }
                value = max(value, sum);
            }
        } else {
            for (IndexType j=1; j<=n; ++j) {
                T sum;
                if (unitDiag) {
                    sum = One;
                    for (IndexType i=j+1; i<=m; ++i) {
                        sum += abs(A(i,j));
                    }
                } else {
                    sum = Zero;
                    for (IndexType i=j; i<=m; ++i) {
                        sum += abs(A(i,j));
                    }
                }
                value = max(value, sum);
            }
        }
    } else if (norm==InfinityNorm) {
//
//      Find normI(A).
//
        if (A.upLo()==Upper) {
            if (A.diag()==Unit) {
                for (IndexType i=1; i<=m; ++i) {
                    work(i) = One;
                }
                for (IndexType j=1; j<=n; ++j) {
                    for (IndexType i=1; i<=min(m,j-1); ++i) {
                        work(i) += abs(A(i,j));
                    }
                }
            } else {
                for (IndexType i=1; i<=m; ++i) {
                    work(i) = Zero;
                }
                for (IndexType j=1; j<=n; ++j) {
                    for (IndexType i=1; i<=min(m,j); ++j) {
                        work(i) += abs(A(i,j));
                    }
                }
            }
        } else {
            if (A.diag()==Unit) {
                for (IndexType i=1; i<=n; ++i) {
                    work(i) = One;
                }
                for (IndexType i=n+1; i<=m; ++i) {
                    work(i) = Zero;
                }
                for (IndexType j=1; j<=n; ++j) {
                    for (IndexType i=j+1; i<=m; ++i) {
                        work(i) += abs(A(i,j));
                    }
                }
            } else {
                for (IndexType i=1; i<=m; ++i) {
                    work(i) = Zero;
                }
                for (IndexType j=1; j<=n; ++j) {
                    for (IndexType i=j; i<=m; ++i) {
                        work(i) += abs(A(i,j));
                    }
                }
            }
        }
        value = Zero;
        for (IndexType i=1; i<=m; ++i) {
            value = max(value, work(i));
        }
    } else if (norm==FrobeniusNorm) {
//
//      Find normF(A).
//
        T  scale, sum;
        if (A.upLo()==Upper) {
            if (A.diag()==Unit) {
                scale = One;
                sum = min(m, n);
                for (IndexType j=2; j<=n; ++j) {
                    const auto rows = _(1,min(m,j-1));
                    const auto Aj = A(rows,j);
                    lassq(Aj, scale, sum);
                }
            } else {
                scale = Zero;
                sum = One;
                for (IndexType j=1; j<=n; ++j) {
                    const auto rows = _(1,min(m,j));
                    const auto Aj = A(rows,j);
                    lassq(Aj, scale, sum);
                }
            }
        } else {
            if (A.diag()==Unit) {
                scale = One,
                sum = min(m,n);
                for (IndexType j=1; j<=n; ++j) {
                    const auto rows = _(min(m,j)+1,m);
                    const auto Aj = A(rows,j);
                    lassq(Aj, scale, sum);
                }
            } else {
                scale = Zero;
                sum = One;
                for (IndexType j=1; j<=n; ++j) {
                    const auto rows = _(min(m+1,j),m);
                    const auto Aj = A(rows,j);
                    lassq(Aj, scale, sum);
                }
            }
        }
        value = scale*sqrt(sum);
    }
    return value;
}

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename VWORK>
typename GeMatrix<MA>::ElementType
lan(Norm norm, const GeMatrix<MA> &A, DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    return cxxlapack::lange<IndexType>(getF77Char(norm),
                                       A.numRows(),
                                       A.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       work.data());
}

template <typename MA, typename VWORK>
typename TrMatrix<MA>::ElementType
lan(Norm norm, const TrMatrix<MA> &A, DenseVector<VWORK> &work)
{
    typedef typename TrMatrix<MA>::IndexType  IndexType;

    return cxxlapack::lantr<IndexType>(getF77Char(norm),
                                       getF77Char(A.upLo()),
                                       getF77Char(A.diag()),
                                       A.numRows(),
                                       A.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       work.data());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//-- lan(ge)
template <typename MA>
typename GeMatrix<MA>::ElementType
lan(Norm norm, const GeMatrix<MA> &A)
{
    ASSERT(norm!=InfinityNorm);

    typedef typename GeMatrix<MA>::ElementType  T;
    DenseVector<Array<T> >  dummy;
    return lan(norm, A, dummy);
}

template <typename MA, typename VWORK>
typename GeMatrix<MA>::ElementType
lan(Norm norm, const GeMatrix<MA> &A, DenseVector<VWORK> &work)
{
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
    T result = lan_generic(norm, A, work);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    if (_work.length()==0) {
        _work.resize(work.length());
    }

    T _result = external::lan(norm, A, _work);

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

//-- lan(tr)
template <typename MA>
typename TrMatrix<MA>::ElementType
lan(Norm norm, const TrMatrix<MA> &A)
{
    ASSERT(norm!=InfinityNorm);

    typedef typename TrMatrix<MA>::ElementType  T;
    DenseVector<Array<T> >  dummy;
    return lan(norm, A, dummy);
}

template <typename MA, typename VWORK>
typename TrMatrix<MA>::ElementType
lan(Norm norm, const TrMatrix<MA> &A, DenseVector<VWORK> &work)
{
    typedef typename TrMatrix<MA>::ElementType T;
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
    T result = lan_generic(norm, A, work);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    if (_work.length()==0) {
        _work.resize(work.length());
    }

    T _result = external::lan(norm, A, _work);

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
lan(Norm norm, const MA &A)
{
    typedef typename MA::ElementType  T;

    CHECKPOINT_ENTER;
    const T result = lan(norm, A);
    CHECKPOINT_LEAVE;

    return result;
}

template <typename MA, typename VWORK>
typename MA::ElementType
lan(Norm norm, const MA &A, VWORK &&work)
{
    typedef typename MA::ElementType  T;

    CHECKPOINT_ENTER;
    const T result = lan(norm, A, work);
    CHECKPOINT_LEAVE;

    return result;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_LAN_TCC
