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

#ifndef FLENS_BLAS_LEVEL1_ROTM_TCC
#define FLENS_BLAS_LEVEL1_ROTM_TCC 1

#include <cxxblas/cxxblas.h>
#include <flens/auxiliary/macros.h>
#include <flens/storage/storage.h>
#include <flens/typedefs.h>

namespace flens { namespace blas {

//-- rotmg
template <typename T, typename VP>
void
rotmg(T &a, T &b, T &c, T &s, DenseVector<VP> &p)
{
#   ifdef HAVE_CXXBLAS_ROTMG
    typedef typename DenseVector<VP>::IndexType  IndexType;

    ASSERT(p.length()==IndexType(5));

    cxxblas::rotmg(p.length(), a, b, c, s, p.data());
#   else
    ASSERT(0);
#   endif
}

//-- forwarding: rotm ----------------------------------------------------------
template <typename VX, typename VY, typename VP>
void
rotm(VX &&x, VY &&y, const VP &p)
{
    rotm(x, y, p);
}

//-- rotm
template <typename VX, typename VY, typename VP>
void
rotm(DenseVector<VX> &x, DenseVector<VY> &y, const DenseVector<VP> &p)
{
#   ifdef HAVE_CXXBLAS_ROTM
    typedef typename DenseVector<VX>::IndexType  IndexType;

    ASSERT(p.length()==IndexType(5));
    ASSERT(x.length()==y.length());

    const IndexType n = x.length();

    cxxblas::rotm(n, x.data(), x.stride(), y.data(), y.stride(), p.data());

#   else
    ASSERT(0);
#   endif
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_ROT_TCC
