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

#ifndef FLENS_BLAS_LEVEL1_AXPY_TCC
#define FLENS_BLAS_LEVEL1_AXPY_TCC 1

#include <cxxblas/cxxblas.h>
#include <flens/aux/macros.h>
#include <flens/blas/debugmacro.h>

namespace flens { namespace blas {

//-- common interface for vectors ----------------------------------------------
template <typename ALPHA, typename VX, typename VY>
void
axpy(const ALPHA &alpha, const Vector<VX> &x, Vector<VY> &y)
{
    axpy(alpha, x.impl(), y.impl());
}

//-- axpy
template <typename ALPHA, typename VX, typename VY>
void
axpy(const ALPHA &alpha, const DenseVector<VX> &x, DenseVector<VY> &y)
{
    FLENS_CLOSURELOG_ADD_ENTRY_AXPY(alpha, x, y);

    if (y.length()==0) {
        y.resize(x.engine(), 0);
    }
    ASSERT(y.length()==x.length());

#   ifdef HAVE_CXXBLAS_AXPY
    cxxblas::axpy(x.length(), alpha,
                  x.data(), x.stride(),
                  y.data(), y.stride());
#   else
    ASSERT(0);
#   endif
    FLENS_CLOSURELOG_END_ENTRY;
}

//-- common interface for matrices ---------------------------------------------
template <typename ALPHA, typename MA, typename MB>
void
axpy(cxxblas::Transpose trans,
     const ALPHA &alpha, const Matrix<MA> &A, Matrix<MB> &B)
{
    axpy(trans, alpha, A.impl(), B.impl());
}

//-- geaxpy
template <typename ALPHA, typename MA, typename MB>
void
axpy(cxxblas::Transpose trans,
     const ALPHA &alpha, const GeMatrix<MA> &A, GeMatrix<MB> &B)
{
    FLENS_CLOSURELOG_ADD_ENTRY_AXPY(alpha, A, B);

    if (B.numRows()*B.numCols()==0) {
        if ((trans==cxxblas::NoTrans) || (trans==cxxblas::Conj)) {
            B.resize(A.numRows(), A.numCols());
        } else {
            B.resize(A.numCols(), A.numRows());
        }
        // fill with zeros!
        B.fill();
    }

#   ifndef NDEBUG
    if ((trans==cxxblas::NoTrans) || (trans==cxxblas::Conj)) {
        ASSERT((A.numRows()==B.numRows()) && (A.numCols()==B.numCols()));
    } else {
        ASSERT((A.numRows()==B.numCols()) && (A.numCols()==B.numRows()));
    }
#   endif

    trans = (MA::order==MB::order)
          ? cxxblas::Transpose(trans ^ cxxblas::NoTrans)
          : cxxblas::Transpose(trans ^ cxxblas::Trans);

#   ifdef HAVE_CXXBLAS_GEAXPY
    cxxblas::geaxpy(MB::order, trans,
                    B.numRows(), B.numCols(), alpha,
                    A.data(), A.leadingDimension(),
                    B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif
    FLENS_CLOSURELOG_END_ENTRY;
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_AXPY_TCC
