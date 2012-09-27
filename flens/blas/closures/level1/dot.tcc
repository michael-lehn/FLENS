/*
 *   Copyright (c) 2007, Michael Lehn
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

#ifndef FLENS_BLAS_CLOSURES_LEVEL1_DOT_TCC
#define FLENS_BLAS_CLOSURES_LEVEL1_DOT_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/closures/closures.h>
#include <flens/blas/level1/level1.h>
#include <flens/blas/level2/level2.h>
#include <flens/blas/level3/level3.h>
#include <flens/typedefs.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

// dot product where x or y is a closure (or unknown vector type)
#ifdef FLENS_DEBUG_CLOSURES

template <typename X, typename Y, typename T>
void
dot(const Vector<X> &x, const Vector<Y> &y, T &result)
{
    FLENS_BLASLOG_BEGIN_DOT(x, y);
//
//  Compute the result of closures x and/or y
//
    typedef typename Vector<X>::Impl   VX;
    typedef typename Vector<Y>::Impl   VY;

    typedef typename Result<VX>::Type  RVX;
    typedef typename Result<VY>::Type  RVY;

    FLENS_BLASLOG_TMP_TRON;

    const RVX &_x = x.impl();
    const RVY &_y = y.impl();

    FLENS_BLASLOG_TMP_TROFF;

//
//  Compute the dot product.  If vectors types of tmpX or tmpY are unknonw
//  we would get an unterminated recursion   So we use a checkpoint guard.
//
    CHECKPOINT_ENTER;
    dot(_x, _y, result);
    CHECKPOINT_LEAVE;

    if (! IsSame<VX,RVX>::value) {
        FLENS_BLASLOG_TMP_REMOVE(_x, x.impl());
    }
    if (! IsSame<VY,RVY>::value) {
        FLENS_BLASLOG_TMP_REMOVE(_y, y.impl());
    }

    FLENS_BLASLOG_END;
}

#else

/*
template <typename X, typename Y, typename T>
void
dot(const Vector<X> &, const Vector<Y> &, T &)
{
//
//  x or y is a closure and we would need a temporary to hold its
//  evaluation.
//
    ASSERT(0);
}
*/

#endif

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_LEVEL1_DOT_TCC
