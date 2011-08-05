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

#ifndef FLENS_BLAS_LEVEL1_SCAL_TCC
#define FLENS_BLAS_LEVEL1_SCAL_TCC 1

#include <cxxblas/cxxblas.h>
#include <flens/aux/macros.h>
#include <flens/storage/storage.h>

namespace flens { namespace blas {

//-- forwarding ----------------------------------------------------------------
template <typename ALPHA, typename VY>
void
scal(const ALPHA &alpha, VY &&y)
{
    scal(alpha, y);
}

//-- common interface for vectors ----------------------------------------------
template <typename ALPHA, typename VY>
void
scal(const ALPHA &alpha, Vector<VY> &y)
{
    scal(alpha, y.impl());
}

//-- common interface for matrices ---------------------------------------------
template <typename ALPHA, typename MB>
void
scal(const ALPHA &alpha, Matrix<MB> &B)
{
    scal(alpha, B.impl());
}

//-- scal
template <typename ALPHA, typename VY>
void
scal(const ALPHA &alpha, DenseVector<VY> &y)
{
#   ifdef HAVE_CXXBLAS_SCAL
    cxxblas::scal(y.length(), alpha, y.data(), y.stride());
#   else
    ASSERT(0);
#   endif
}

//-- gescal
template <typename ALPHA, typename MB>
void
scal(const ALPHA &alpha, GeMatrix<MB> &B)
{
#   ifdef HAVE_CXXBLAS_GESCAL
    cxxblas::gescal(StorageInfo<MB>::Order,
                    B.numRows(), B.numCols(),
                    alpha, B.engine().data(), B.engine().leadingDimension());
#   else
    ASSERT(0);
#   endif
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_SCAL_TCC
