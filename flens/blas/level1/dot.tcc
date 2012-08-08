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

#ifndef FLENS_BLAS_LEVEL1_DOT_TCC
#define FLENS_BLAS_LEVEL1_DOT_TCC 1

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

template <typename X, typename Y, typename T>
void
dot(const DenseVector<X> &x, const DenseVector<Y> &y, T &result)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_DOT(x, y);

    ASSERT(x.length()==y.length());

#   ifdef HAVE_CXXBLAS_DOT
    cxxblas::dot(x.length(),
                 x.data(), x.stride(),
                 y.data(), y.stride(), result);
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

template <typename X, typename Y, typename T>
void
dotu(const DenseVector<X> &x, const DenseVector<Y> &y, T &result)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_DOTU(x, y);

    ASSERT(x.length()==y.length());

#   ifdef HAVE_CXXBLAS_DOTU
    cxxblas::dotu(x.length(),
                  x.data(), x.stride(),
                  y.data(), y.stride(), result);
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

template <typename X, typename Y>
typename RestrictTo<IsSame<typename X::ElementType,
                           typename Y::ElementType>::value,
                           typename X::ElementType>::Type
dot(const DenseVector<X> &x, const DenseVector<Y> &y)
{
    typename X::ElementType val;
    dot(x, y, val);
    return val;
}

template <typename X, typename Y>
typename RestrictTo<IsSame<typename X::ElementType,
                           typename Y::ElementType>::value,
                           typename X::ElementType>::Type
dotu(const DenseVector<X> &x, const DenseVector<Y> &y)
{
    typename X::ElementType val;
    dotu(x, y, val);
    return val;
}


} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_DOT_TCC
