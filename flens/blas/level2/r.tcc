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

#ifndef FLENS_BLAS_LEVEL2_R_TCC
#define FLENS_BLAS_LEVEL2_R_TCC 1

#include <cxxblas/cxxblas.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//-- forwarding ----------------------------------------------------------------
template <typename ALPHA, typename VX, typename VY, typename MA>
void
r(const ALPHA &alpha, const VX &x, const VY &y, MA &&A)
{
    CHECKPOINT_ENTER;
    r(alpha, x, y, A);
    CHECKPOINT_LEAVE;
}

template <typename ALPHA, typename VX, typename VY, typename MA>
void
ru(const ALPHA &alpha, const VX &x, const VY &y, MA &&A)
{
    CHECKPOINT_ENTER;
    ru(alpha, x, y, A);
    CHECKPOINT_LEAVE;
}

template <typename ALPHA, typename VX, typename VY, typename MA>
void
rc(const ALPHA &alpha, const VX &x, const VY &y, MA &&A)
{
    CHECKPOINT_ENTER;
    rc(alpha, x, y, A);
    CHECKPOINT_LEAVE;
}

//-- GeMatrix, DenseVector -----------------------------------------------------

//-- ger
template <typename ALPHA, typename VX, typename VY, typename MA>
void
r(const ALPHA &alpha,
  const DenseVector<VX> &x, const DenseVector<VY> &y,
  GeMatrix<MA> &A)
{
    if ((x.length()!=A.numRows()) || (y.length()!=A.numCols())) {
        A.resize(x.length(), y.length(), x.firstIndex(), y.firstIndex());
    }
    cxxblas::ger(A.order(),
                 A.numRows(), A.numCols(),
                 alpha,
                 x.data(), x.stride(),
                 y.data(), y.stride(),
                 A.data(), A.leadingDimension());
}

//-- geru
template <typename ALPHA, typename VX, typename VY, typename MA>
void
ru(const ALPHA &alpha,
   const DenseVector<VX> &x, const DenseVector<VY> &y,
   GeMatrix<MA> &A)
{
    r(alpha, x, y, A);
}

//-- gerc
template <typename ALPHA, typename VX, typename VY, typename MA>
void
rc(const ALPHA &alpha,
   const DenseVector<VX> &x, const DenseVector<VY> &y,
   GeMatrix<MA> &A)
{
    if ((x.length()!=A.numRows()) || (y.length()!=A.numCols())) {
        A.resize(x.length(), y.length(), x.firstIndex(), y.firstIndex());
    }
    cxxblas::gerc(A.order(),
                  A.numRows(), A.numCols(),
                  alpha,
                  x.data(), x.stride(),
                  y.data(), y.stride(),
                  A.data(), A.leadingDimension());
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL2_R_TCC
