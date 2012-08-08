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

#include <cxxblas/cxxblas.h>
#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/closures/closures.h>
#include <flens/blas/level1/level1.h>
#include <flens/typedefs.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif
namespace flens { namespace blas {

//-- rotg
template <typename T>
void
rotg(T &a, T &b, T &c, T &s)
{
#   ifdef HAVE_CXXBLAS_ROTG
    rotg(a, b, c, s);
#   else
    ASSERT(0);
#   endif
}


//-- forwarding: rot -----------------------------------------------------------
template <typename VX, typename VY, typename T>
void
rot(VX &&x, VY &&y, const T &c, const T &s)
{
    rot(x, y, c, s);
}

//-- rot
template <typename VX, typename VY, typename T>
void
rot(DenseVector<VX> &x, DenseVector<VY> &y, const T &c, const T &s)
{
    ASSERT(x.length()==y.length());
#   ifdef HAVE_CXXBLAS_ROT
    cxxblas::rot(x.length(), x.data(), x.stride(), y.data(), y.stride(), c, s);
#   else
    ASSERT(0);
#   endif
}


} } // namespace blas, flens
