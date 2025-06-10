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

//-- rotmg
template <typename T, typename VP>
void
rotmg(T &d1, T &d2, T &b1, T &b2, DenseVector<VP> &p)
{
#   ifdef HAVE_CXXBLAS_ROTMG

#   ifndef NDEBUG
    typedef typename DenseVector<VP>::IndexType  IndexType;
#   endif

    ASSERT(p.length()==IndexType(5));

    cxxblas::rotmg(d1, d2, b1, b2, p.data());
#   else
    ASSERT(0);
    (void)(d1);
    (void)(d2);
    (void)(b1);
    (void)(b2);
    (void)(p);
#   endif
}

//-- rotm
template <typename VX, typename VY, typename VP>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
rotm(VX &&x, VY &&y, const DenseVector<VP> &p)
{
#   ifdef HAVE_CXXBLAS_ROTM
    typedef typename RemoveRef<VX>::Type    VectorX;
    typedef typename VectorX::IndexType     IndexType;

    ASSERT(p.length()==IndexType(5));
    ASSERT(x.length()==y.length());

    const IndexType n = x.length();

    cxxblas::rotm(n, x.data(), x.stride(), y.data(), y.stride(), p.data());

#   else
    ASSERT(0);
    (void)(x);
    (void)(y);
    (void)(p);
#   endif
}

//-- rotmg
template <typename T, typename VP>
void
rotmg(T &d1, T &d2, T &b1, T &b2, TinyVector<VP> &p)
{
#   ifdef HAVE_CXXBLAS_ROTMG

#   ifndef NDEBUG
    typedef typename TinyVector<VP>::IndexType  IndexType;
#   endif

    const int length = VP::length;
    ASSERT(length == IndexType(5));

    cxxblas::rotmg(d1, d2, b1, b2, p.data());
#   else
    ASSERT(0);
    (void)(d1);
    (void)(d2);
    (void)(b1);
    (void)(b2);
    (void)(p);
#   endif
}

//-- rotm
template <typename VX, typename VY, typename VP>
typename RestrictTo<IsTinyVector<VX>::value
                 && IsTinyVector<VY>::value,
         void>::Type
rotm(VX &&x, VY &&y, const TinyVector<VP> &p)
{
#   ifdef HAVE_CXXBLAS_ROTM
    typedef typename std::remove_reference<VX>::type VectorX;
    typedef typename VectorX::ElementType          TX;
    typedef typename std::remove_reference<VY>::type VectorY;
    typedef typename VectorY::ElementType          TY;

    const int lengthX = VectorX::Engine::length;
    const int strideX = VectorX::Engine::stride;
    const int lengthY = VectorY::Engine::length;
    const int strideY = VectorY::Engine::stride;

    const int length = VP::length;

    typedef typename VP::IndexType IndexType;

    ASSERT(length == IndexType(5));
    ASSERT(lengthX == lengthY);

    cxxblas::rotm(lengthX, x.data(), strideX, y.data(), strideY, p.data());

#   else
    ASSERT(0);
    (void)(x);
    (void)(y);
    (void)(p);
#   endif
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_ROT_TCC
