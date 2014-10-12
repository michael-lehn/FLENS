/*
 *   Copyright (c) 2014, Michael Lehn
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
       SUBROUTINE ZLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LASR_TCC
#define FLENS_LAPACK_LA_LASR_TCC 1

#include <cxxstd/cmath.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename VC, typename VS, typename MA>
void
lasr_impl(Side              side,
          LASR::Pivot       pivot,
          LASR::Direction   direction,
          DenseVector<VC>   &c,
          DenseVector<VS>   &s,
          GeMatrix<MA>      &A)
{
    using std::max;

    typedef typename MA::IndexType    IndexType;
    typedef typename MA::ElementType  T;

    const T Zero(0), One(1);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

//
//     Quick return if possible
//
    if (m==0 || n==0) {
        return;
    }

    if (side==Left) {
//
//      Form  P * A
//
        if (pivot==LASR::VariablePivot) {
            if (direction==LASR::Forward) {
                for (IndexType j=1; j<=m-1; ++j) {
                    const T cTmp = c(j);
                    const T sTmp = s(j);
                    if (cTmp!=One || sTmp!=Zero) {
                        for (IndexType i=1; i<=n; ++i) {
                            const T tmp = A(j+1,i);
                            A(j+1, i) = cTmp*tmp - sTmp*A(j, i);
                            A(j,   i) = sTmp*tmp + cTmp*A(j, i);
                        }
                    }
                }
            } else if (direction==LASR::Backward) {
                for (IndexType j=m-1; j>=1; --j) {
                    const T cTmp = c(j);
                    const T sTmp = s(j);
                    if (cTmp!=One || sTmp!=Zero) {
                        for (IndexType i=1; i<=n; ++i) {
                            const T tmp = A(j+1, i);
                            A(j+1, i) = cTmp*tmp - sTmp*A(j, i);
                            A(j,   i) = sTmp*tmp + cTmp*A(j, i);
                        }
                    }
                }
            }
        } else if (pivot==LASR::TopPivot) {
            if (direction==LASR::Forward) {
                for (IndexType j=2; j<=m; ++j) {
                    const T cTmp = c(j-1);
                    const T sTmp = s(j-1);
                    if (cTmp!=One || sTmp!=Zero) {
                        for (IndexType i=1; i<=n; ++i) {
                            const T tmp = A(j, i);
                            A(j, i) = cTmp*tmp - sTmp*A(1, i);
                            A(1, i) = sTmp*tmp + cTmp*A(1, i);
                        }
                    }
                }
            } else if (direction==LASR::Backward) {
                for (IndexType j=m; j>=2; --j) {
                    const T cTmp = c(j-1);
                    const T sTmp = s(j-1);
                    if (cTmp!=One || sTmp!=Zero) {
                        for (IndexType i=1; i<=n; ++i) {
                            const T tmp = A(j, i);
                            A(j, i) = cTmp*tmp - sTmp*A(1, i);
                            A(1, i) = sTmp*tmp + cTmp*A(1, i);
                        }
                    }
                }
            }
        } else if (pivot==LASR::BottomPivot) {
            if (direction==LASR::Forward) {
                for (IndexType j=1; j<=m-1; ++j) {
                    const T cTmp = c(j);
                    const T sTmp = s(j);
                    if (cTmp!=One || sTmp!=Zero) {
                        for (IndexType i=1; i<=n; ++i) {
                            const T tmp = A(j, i);
                            A(j, i) = sTmp*A(m, i) + cTmp*tmp;
                            A(m, i) = cTmp*A(m, i) - sTmp*tmp;
                        }
                    }
                }
            } else if (direction==LASR::Backward) {
                for (IndexType j=m-1; j>=1; --j) {
                    const T cTmp = c(j);
                    const T sTmp = s(j);
                    if (cTmp!=One || sTmp!=Zero) {
                        for (IndexType i=1; i<=n; ++i) {
                            const T tmp = A(j, i);
                            A(j, i) = sTmp*A(m, i) + cTmp*tmp;
                            A(m, i) = cTmp*A(m, i) - sTmp*tmp;
                        }
                    }
                }
            }
        }
    } else if (side==Right) {
//
//      Form A * P**T
//
        if (pivot==LASR::VariablePivot) {
            if (direction==LASR::Forward) {
                for (IndexType j=1; j<=n-1; ++j) {
                    const T cTmp = c(j);
                    const T sTmp = s(j);
                    if (cTmp!=One || sTmp!=Zero) {
                        for (IndexType i=1; i<=m; ++i) {
                            const T tmp = A(i, j+1);
                            A(i, j+1) = cTmp*tmp - sTmp*A(i, j);
                            A(i, j  ) = sTmp*tmp + cTmp*A(i, j);
                        }
                    }
                }
            } else if (direction==LASR::Backward) {
                for (IndexType j=n-1; j>=1; --j) {
                    const T cTmp = c(j);
                    const T sTmp = s(j);
                    if (cTmp!=One || sTmp!=Zero) {
                        for (IndexType i=1; i<=m; ++i) {
                            const T tmp = A( i, j+1);
                            A(i, j+1) = cTmp*tmp - sTmp*A(i, j);
                            A(i, j  ) = sTmp*tmp + cTmp*A(i, j);
                        }
                    }
                }
            }
        } else if (pivot==LASR::TopPivot) {
            if (direction==LASR::Forward) {
                for (IndexType j=2; j<=n; ++j) {
                    const T cTmp = c(j-1);
                    const T sTmp = s(j-1);
                    if (cTmp!=One || sTmp!=Zero) {
                        for (IndexType i=1; i<=m; ++i) {
                            const T tmp = A(i, j);
                            A(i, j) = cTmp*tmp - sTmp*A(i, 1);
                            A(i, 1) = sTmp*tmp + cTmp*A(i, 1);
                        }
                    }
                }
            } else if (direction==LASR::Backward) {
                for (IndexType j=n; j>=2; --j) {
                    const T cTmp = c(j-1);
                    const T sTmp = s(j-1);
                    if (cTmp!=One || sTmp!=Zero) {
                        for (IndexType i=1; i<=m; ++i) {
                            const T tmp = A(i, j);
                            A(i, j) = cTmp*tmp - sTmp*A(i, 1);
                            A(i, 1) = sTmp*tmp + cTmp*A(i, 1);
                        }
                    }
                }
            }
        } else if (pivot==LASR::BottomPivot) {
            if (direction==LASR::Forward) {
                for (IndexType j=1; j<=n-1; ++j) {
                    const T cTmp = c(j);
                    const T sTmp = s(j);
                    if (cTmp!=One || sTmp!=Zero) {
                        for (IndexType i=1; i<=m; ++i) {
                            const T tmp = A(i, j);
                            A(i, j) = sTmp*A(i, n) + cTmp*tmp;
                            A(i, n) = cTmp*A(i, n) - sTmp*tmp;
                        }
                    }
                }
            } else if (direction==LASR::Backward) {
                for (IndexType j=n-1; j>=1; --j) {
                    const T cTmp = c(j);
                    const T sTmp = s(j);
                    if (cTmp!=One || sTmp!=Zero) {
                        for (IndexType i=1; i<=m; ++i) {
                            const T tmp = A(i, j);
                            A(i, j) = sTmp*A(i, n) + cTmp*tmp;
                            A(i, n) = cTmp*A(i, n) - sTmp*tmp;
                        }
                    }
                }
            }
        }
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename VC, typename VS, typename MA>
void
lasr_impl(Side              side,
          LASR::Pivot       pivot,
          LASR::Direction   direction,
          DenseVector<VC>   &c,
          DenseVector<VS>   &s,
          GeMatrix<MA>      &A)
{
    cxxlapack::lasr(getF77Char(side),
                    getF77Char(pivot),
                    getF77Char(direction),
                    A.numRows(),
                    A.numCols(),
                    c.data(),
                    s.data(),
                    A.data(),
                    A.leadingDimension());
}

} // namespace external

#endif

//== public interface ==========================================================

template <typename VC, typename VS, typename MA>
typename RestrictTo<IsRealDenseVector<VC>::value
                 && IsRealDenseVector<VS>::value
                 && IsComplexGeMatrix<MA>::value,
         void>::Type
lasr(Side             side,
     LASR::Pivot      pivot,
     LASR::Direction  direction,
     VC               &&c,
     VS               &&s,
     MA               &&A)
{
    LAPACK_DEBUG_OUT("lasr");
//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<VC>::Type    VectorC;
    typedef typename RemoveRef<VS>::Type    VectorS;
#   endif

#   if defined CHECK_CXXLAPACK || !defined NDEBUG
    typedef typename RemoveRef<MA>::Type    MatrixA;
#   endif
//
//  Test the input parameters
//
#   ifndef NDEBUG

    typedef typename MatrixA::IndexType  IndexType;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    ASSERT(c.firstIndex()==1);
    ASSERT(s.firstIndex()==1);
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    if (side==Left) {
        ASSERT(c.length()==m-1);
        ASSERT(s.length()==m-1);
    } else {
        ASSERT(c.length()==n-1);
        ASSERT(s.length()==n-1);
    }
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename VectorC::NoView      c_org = c;
    typename VectorS::NoView      s_org = s;
    typename MatrixA::NoView      A_org = A;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::lasr_impl(side, pivot, direction, c, s, A);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename VectorC::NoView      c_generic = c;
    typename VectorS::NoView      s_generic = s;
    typename MatrixA::NoView      A_generic = A;

    c = c_org;
    s = s_org;
    A = A_org;

    external::lasr_impl(side, pivot, direction, c, s, A);

    bool failed = false;
    if (! isIdentical(c_generic, c, "c_generic", "c")) {
        std::cerr << "CXXLAPACK: c_generic = " << c_generic << std::endl;
        std::cerr << "F77LAPACK: c = " << c << std::endl;
    }
    if (! isIdentical(s_generic, s, "s_generic", "s")) {
        std::cerr << "CXXLAPACK: s_generic = " << s_generic << std::endl;
        std::cerr << "F77LAPACK: s = " << s << std::endl;
    }
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
    }
    if (failed) {
        ASSERT(0);
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LASR_TCC
