/*
 *   Copyright (c) 2009, Michael Lehn
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

#ifndef FLENS_BLAS_LEVEL1_COPY_TCC
#define FLENS_BLAS_LEVEL1_COPY_TCC 1

#include <flens/aux/macros.h>
#include <flens/blas/blas.h>
#include <flens/blas/debugmacro.h>
#include <flens/typedefs.h>

namespace flens { namespace blas {

namespace impl {

//== specialized interfaces for vectors ========================================

//-- copy
template <typename VX, typename VY>
void
copy(const DenseVector<VX> &x, DenseVector<VY> &y)
{
    FLENS_CLOSURELOG_ADD_ENTRY_COPY(x, y);
    if (y.length()!=x.length()) {
        y.resize(x);
    }
    y.changeIndexBase(x.firstIndex());

#   ifdef HAVE_CXXBLAS_COPY
    cxxblas::copy(x.length(), x.data(), x.stride(), y.data(), y.stride());
#   else
    ASSERT(0);
#   endif

    FLENS_CLOSURELOG_END_ENTRY
}

//== specialized interfaces for matrices =======================================

//-- gecopy
template <typename MA, typename MB>
void
copy(Transpose trans, const GeMatrix<MA> &A, GeMatrix<MB> &B)
{
    FLENS_CLOSURELOG_ADD_ENTRY_COPY(A, B);

    if (trans==NoTrans) {
        if ((A.numRows()!=B.numRows()) || (A.numCols()!=B.numCols())) {
            B.resize(A);
        }
    } else {
        if ((A.numRows()!=B.numCols())  || (A.numCols()!=B.numRows())) {
            B.resize(A.numCols(), A.numRows(),
                     A.firstCol(), A.firstRow());
        }
    }
    trans = (MA::order==MB::order)
          ? Transpose(trans ^ NoTrans)
          : Transpose(trans ^ Trans);

#   ifdef HAVE_CXXBLAS_GECOPY
    cxxblas::gecopy(MB::order, trans,
                    B.numRows(), B.numCols(),
                    A.data(), A.leadingDimension(),
                    B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif
    FLENS_CLOSURELOG_END_ENTRY
}

//-- trcopy
template <typename MA, typename MB>
void
copy(Transpose trans, const TrMatrix<MA> &A, TrMatrix<MB> &B)
{
    FLENS_CLOSURELOG_ADD_ENTRY_COPY(A, B);

    // TODO: change default resize policy such that resize is done:
    //          1) if possible, e.g. no resize in an axpy
    //          2) rhs-size is zero
    if (B.dim()!=A.dim() && B.dim()==0) {
        B.resize(A);
        B.upLo() = A.upLo();
        B.diag() = A.diag();
    }

#   ifndef NDEBUG
    if (trans==NoTrans) {
        ASSERT(A.upLo()==B.upLo());
    } else {
        ASSERT(A.upLo()!=B.upLo());
    }
#   endif

    // TODO: make this assertion unnecessary
    ASSERT(A.order()==B.order());
    ASSERT(A.diag()==B.diag());

#   ifdef HAVE_CXXBLAS_TRCOPY
    cxxblas::trcopy(B.order(), B.upLo(), trans, B.diag(), B.dim(),
                    A.data(), A.leadingDimension(),
                    B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif
    FLENS_CLOSURELOG_END_ENTRY
}

//-- sycopy
template <typename MA, typename MB>
void
copy(const SyMatrix<MA> &A, SyMatrix<MB> &B)
{
    FLENS_CLOSURELOG_ADD_ENTRY_COPY(A, B);

    // TODO: change default resize policy such that resize is done:
    //          1) if possible, e.g. no resize in an axpy
    //          2) rhs-size is zero
    if (B.dim()!=A.dim() && B.dim()==0) {
        B.resize(A);
        B.upLo() = A.upLo();
    }

    // TODO: make this assertion unnecessary
    ASSERT(A.order()==B.order());

#   ifdef HAVE_CXXBLAS_GECOPY
    cxxblas::trcopy(B.order(), B.upLo(), NoTrans, NonUnit, B.dim(),
                    A.data(), A.leadingDimension(),
                    B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif
    FLENS_CLOSURELOG_END_ENTRY
}

//== common interface for vectors ==============================================
template <typename VX, typename VY>
void
copy(const Vector<VX> &x, Vector<VY> &y)
{
    CHECKPOINT_ENTER;
    copy(x.impl(), y.impl());
    CHECKPOINT_LEAVE;
}

//== common interface for matrices =============================================
template <typename MA, typename MB>
void
copy(Transpose trans, const Matrix<MA> &A, Matrix<MB> &B)
{
    CHECKPOINT_ENTER;
    copy(trans, A.impl(), B.impl());
    CHECKPOINT_LEAVE;
}

template <typename MA, typename MB>
void
copy(Transpose trans, const SymmetricMatrix<MA> &A, SymmetricMatrix<MB> &B)
{
    //TODO: Lehn: this should become a warning instead of an error
    ASSERT(trans==NoTrans);

    CHECKPOINT_ENTER;
    copy(A.impl(), B.impl());
    CHECKPOINT_LEAVE;
}

} // namespace impl

//== entry points ==============================================================
template <typename X, typename Y>
void
copy(const X &x, Y &&y)
{
    CHECKPOINT_ENTER;
    impl::copy(x.impl(), y.impl());
    CHECKPOINT_LEAVE;
}

template <typename MA, typename MB>
void
copy(Transpose trans, const MA &A, MB &&B)
{
    CHECKPOINT_ENTER;
    impl::copy(trans, A.impl(), B.impl());
    CHECKPOINT_LEAVE;
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_COPY_TCC
