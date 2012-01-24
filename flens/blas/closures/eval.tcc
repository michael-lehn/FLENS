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

#ifndef FLENS_BLAS_CLOSURES_EVAL_TCC
#define FLENS_BLAS_CLOSURES_EVAL_TCC 1

#include <flens/blas/closures/debugclosure.h>
#include <flens/blas/closures/prunematrixclosure.h>
#include <flens/blas/closures/prunevectorclosure.h>
#include <flens/blas/closures/result.h>
#include <flens/blas/debugmacro.h>
#include <flens/blas/level1/level1.h>
#include <flens/blas/level2/level2.h>
#include <flens/typedefs.h>


//-- entry points---------------------------------------------------------------

namespace flens {

//-- vector closures
template <typename VX, typename VY>
void
assign(const Vector<VX> &x, Vector<VY> &y)
{
    CHECKPOINT_ENTER;
    FLENS_CLOSURELOG_BEGIN_ASSIGNMENT(x, y);

    blas::copy(x.impl(), y.impl());

    FLENS_CLOSURELOG_END;
    CHECKPOINT_LEAVE;
}

} // namespace flens


namespace flens { namespace blas {

//-- vector closures -----------------------------------------------------------
//
// y = x1 + x2
//
template <typename VL, typename VR, typename VY>
void
copy(const VectorClosure<OpAdd, VL, VR> &x, Vector<VY> &y)
{
    CHECKPOINT_ENTER;
    FLENS_CLOSURELOG_BEGIN_COPY(x, y);

    typedef typename VY::Impl::ElementType T;
    const T  One(1);
//
//  y = x1
//
    blas::copy(x.left(), y.impl());
//
//  y += x2
//
    blas::axpy(One, x.right(), y.impl());

    FLENS_CLOSURELOG_END;
    CHECKPOINT_LEAVE;
}

//
// y = x1 - x2
//
template <typename VL, typename VR, typename VY>
void
copy(const VectorClosure<OpSub, VL, VR> &x, Vector<VY> &y)
{
    CHECKPOINT_ENTER;
    FLENS_CLOSURELOG_BEGIN_COPY(x, y);

    typedef typename VY::Impl::ElementType T;
    const T  MinusOne(-1);
//
//  y = x1
//
    blas::copy(x.left(), y.impl());
//
//  y -= x2
//
    blas::axpy(MinusOne, x.right(), y.impl());

    FLENS_CLOSURELOG_END;
    CHECKPOINT_LEAVE;
}

//-- axpy
// y += alpha*(x1+x2)
//
//
template <typename ALPHA, typename VL, typename VR, typename VY>
void
axpy(const ALPHA &alpha,
     const VectorClosure<OpAdd, VL, VR> &x, Vector<VY> &y)
{
    typedef VectorClosure<OpAdd, VL, VR>  VC;

#   ifndef FLENS_DEBUG_CLOSURES
//
//  Note that y += x1+x2  or  y -= x1+x2 are equivalent to y = y + (x1+x2) or
//  y = y - (x1+x2) respectively.  Hence, the result of x1+x2 has to be
//  computed *before* we can update y.  This means that a temporary vectors is
//  needed to hold the result of x1+x2.  We only allow this if the
//  FLENS_DEBUG_CLOSURES macro is defined such that a user has a chance to
//  optimize his/her expressions.
//
    ASSERT(0);

#   else
    CHECKPOINT_ENTER;
    FLENS_CLOSURELOG_BEGIN_AXPY(alpha, x, y);

//
//  Compute the result of closure x = (x.left()+x.right()) first and store
//  it in tmp.
//
    typename Result<VC>::Type tmp;
    FLENS_CLOSURELOG_ADD_TMP(tmp, x);

    tmp = x;
//
//  Update y with tmp, i.e. compute y = y + alpha*tmp
//
    axpy(alpha, tmp, y.impl());

    FLENS_CLOSURELOG_END;
    CHECKPOINT_LEAVE;
#   endif
}

//
// y += x1 - x2  or  y -= x1 - x2  so alpha = 1 or alpha = -1
// otherwise this would be y += alpha*(x1-x2) and requires a
// temporary for x1-x2.
//
template <typename ALPHA, typename VL, typename VR, typename VY>
void
axpy(const ALPHA &alpha,
     const VectorClosure<OpSub, VL, VR> &x, Vector<VY> &y)
{
    typedef VectorClosure<OpSub, VL, VR>  VC;

#   ifndef FLENS_DEBUG_CLOSURES
//
//  Note that y += x1-x2  or  y -= x1-x2 are equivalent to y = y + (x1-x2) or
//  y = y - (x1-x2) respectively.  Hence, the result of x1-x2 has to be
//  computed *before* we can update y.  This means that a temporary vectors is
//  needed to hold the result of x1-x2.  We only allow this if the
//  FLENS_DEBUG_CLOSURES macro is defined such that a user has a chance to
//  optimize his/her expressions.
//
    ASSERT(0);

#   else
    CHECKPOINT_ENTER;
    FLENS_CLOSURELOG_BEGIN_AXPY(alpha, x, y);

//
//  Compute the result of closure x = (x.left()-x.right()) first and store
//  it in tmp.
//
    typename Result<VC>::Type tmp;
    FLENS_CLOSURELOG_ADD_TMP(tmp, x);

    tmp = x;
//
//  Update y with tmp, i.e. compute y = y + alpha*tmp
//
    axpy(alpha, tmp, y.impl());

    FLENS_CLOSURELOG_END;
    CHECKPOINT_LEAVE;
#   endif
}

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_EVAL_TCC
