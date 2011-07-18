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

#ifndef FLENS_SCALAROPERATIONS_H
#define FLENS_SCALAROPERATIONS_H 1

#include <flens/matvec.h>
#include <flens/operationtypes.h>
#include <flens/scalarclosures.h>
#include <flens/traits.h>

namespace flens {

//== auxiliary traits for ScalarClosures =======================================

//-- trait for the compuation of scalar operations -----------------------------
template <typename Op>
struct Operation
{
};

//-- trait for ScalarOperation -------------------------------------------------
template <typename Op, typename X, typename Y>
    struct _SO_;

//-- trait for ScalarUnaryOperation --------------------------------------------

template <typename Op, typename X>
    struct _SUO_;

//==============================================================================

//-- x+y -----------------------------------------------------------------------
template <typename X, typename Y>
    const typename _SO_<OpAdd, X, Y>::Type
    operator+(const ScalarClosure<X> &x, const ScalarClosure<Y> &y);

template <typename X, typename T>
    const typename _SO_<OpAdd, X, Scalar<T> >::Type
    operator+(const ScalarClosure<X> &x, const T &a);

template <typename T, typename X>
    const typename _SO_<OpAdd, Scalar<T>, X>::Type
    operator+(const T &a, const ScalarClosure<X> &x);

//-- x-y -----------------------------------------------------------------------
template <typename X, typename Y>
    const typename _SO_<OpMinus, X, Y>::Type
    operator-(const ScalarClosure<X> &x, const ScalarClosure<Y> &y);

template <typename X, typename T>
    const typename _SO_<OpMinus, X, Scalar<T> >::Type
    operator-(const ScalarClosure<X> &x, const T &a);

template <typename T, typename X>
    const typename _SO_<OpMinus, Scalar<T>, X>::Type
    operator-(const T &a, const ScalarClosure<X> &x);

//-- x*y -----------------------------------------------------------------------
template <typename X, typename Y>
    const typename _SO_<OpMult, X, Y>::Type
    operator*(const ScalarClosure<X> &x, const ScalarClosure<Y> &y);

template <typename X, typename T>
    const typename _SO_<OpMult, X, Scalar<T> >::Type
    operator*(const ScalarClosure<X> &x, const T &a);

template <typename T, typename X>
    const typename _SO_<OpMult, Scalar<T>, X>::Type
    operator*(const T &a, const ScalarClosure<X> &x);

//-- x/y -----------------------------------------------------------------------
template <typename X, typename Y>
    const typename _SO_<OpDiv, X, Y>::Type
    operator/(const ScalarClosure<X> &x, const ScalarClosure<Y> &y);

template <typename X, typename T>
    const typename _SO_<OpDiv, X, Scalar<T> >::Type
    operator/(const ScalarClosure<X> &x, const T &a);

template <typename T, typename X>
    const typename _SO_<OpDiv, Scalar<T>, X>::Type
    operator/(const T &a, const ScalarClosure<X> &x);

//-- x%y -----------------------------------------------------------------------
template <typename X, typename Y>
    const typename _SO_<OpMod, X, Y>::Type
    operator%(const ScalarClosure<X> &x, const ScalarClosure<Y> &y);

template <typename X, typename T>
    const typename _SO_<OpMod, X, Scalar<T> >::Type
    operator%(const ScalarClosure<X> &x, const T &a);

template <typename T, typename X>
    const typename _SO_<OpMod, Scalar<T>, X>::Type
    operator%(const T &a, const ScalarClosure<X> &x);

//-- max(x,y) ------------------------------------------------------------------
template <typename X, typename Y>
    const typename _SO_<OpMax, X, Y>::Type
    max(const ScalarClosure<X> &x, const ScalarClosure<Y> &y);

template <typename X>
    const typename _SO_<OpMax, X, X>::Type
    max(const ScalarClosure<X> &x, const ScalarClosure<X> &y);

template <typename X, typename T>
    const typename _SO_<OpMax, X, Scalar<T> >::Type
    max(const ScalarClosure<X> &x, const T &a);

template <typename T, typename X>
    const typename _SO_<OpMax, Scalar<T>, X>::Type
    max(const T &a, const ScalarClosure<X> &x);

//-- min(x,y) ------------------------------------------------------------------
template <typename X, typename Y>
    const typename _SO_<OpMin, X, Y>::Type
    min(const ScalarClosure<X> &x, const ScalarClosure<Y> &y);

template <typename X>
    const typename _SO_<OpMin, X, X>::Type
    min(const ScalarClosure<X> &x, const ScalarClosure<X> &y);

template <typename X, typename T>
    const typename _SO_<OpMin, X, Scalar<T> >::Type
    min(const ScalarClosure<X> &x, const T &a);

template <typename T, typename X>
    const typename _SO_<OpMin, Scalar<T>, X>::Type
    min(const T &a, const ScalarClosure<X> &x);

//-- pow(x,y) ------------------------------------------------------------------
template <typename X, typename Y>
    const typename _SO_<OpPow, X, Y>::Type
    pow(const ScalarClosure<X> &x, const ScalarClosure<Y> &y);

template <typename X>
    const typename _SO_<OpPow, X, X>::Type
    pow(const ScalarClosure<X> &x, const ScalarClosure<X> &y);

template <typename X, typename T>
    const typename _SO_<OpPow, X, Scalar<T> >::Type
    pow(const ScalarClosure<X> &x, const T &a);

template <typename T, typename X>
    const typename _SO_<OpPow, Scalar<T>, X>::Type
    pow(const T &a, const ScalarClosure<X> &x);

//-- -x ------------------------------------------------------------------------
template <typename X>
    const typename _SUO_<OpMinus, X>::Type
    operator-(const ScalarClosure<X> &x);

//-- conjugate(x), conj(x) -----------------------------------------------------
template <typename X>
    const typename _SUO_<OpConj, X>::Type
    conjugate(const ScalarClosure<X> &x);

template <typename X>
    const typename _SUO_<OpConj, X>::Type
    conj(const ScalarClosure<X> &x);

//-- ceil(x) -------------------------------------------------------------------
template <typename X>
    const typename _SUO_<OpCeil, X>::Type
    ceil(const ScalarClosure<X> &x);

//-- floor(x) ------------------------------------------------------------------
template <typename X>
    const typename _SUO_<OpFloor, X>::Type
    floor(const ScalarClosure<X> &x);

//-- abs(x) --------------------------------------------------------------------
template <typename X, typename Y>
    const typename _SUO_<OpAbs, X>::Type
    abs(const ScalarClosure<X> &x);

//-- sqrt(x) -------------------------------------------------------------------
template <typename X>
    const typename _SUO_<OpSqrt, X>::Type
    sqrt(const ScalarClosure<X> &x);

//-- exp(x) --------------------------------------------------------------------
template <typename X>
    const typename _SUO_<OpExp, X>::Type
    exp(const ScalarClosure<X> &x);

//-- log(x) --------------------------------------------------------------------
template <typename X>
    const typename _SUO_<OpLog, X>::Type
    log(const ScalarClosure<X> &x);

//-- cos(x) --------------------------------------------------------------------
template <typename X>
    const typename _SUO_<OpCos, X>::Type
    cos(const ScalarClosure<X> &x);

//-- sin(x) --------------------------------------------------------------------
template <typename X>
    const typename _SUO_<OpSin, X>::Type
    sin(const ScalarClosure<X> &x);

//-- tan(x) --------------------------------------------------------------------
template <typename X>
    const typename _SUO_<OpTan, X>::Type
    tan(const ScalarClosure<X> &x);

//-- cosh(x) -------------------------------------------------------------------
template <typename X>
    const typename _SUO_<OpCosh, X>::Type
    cosh(const ScalarClosure<X> &x);

//-- sinh(x) -------------------------------------------------------------------
template <typename X>
    const typename _SUO_<OpSinh, X>::Type
    sinh(const ScalarClosure<X> &x);

//-- tanh(x) -------------------------------------------------------------------
template <typename X>
    const typename _SUO_<OpTanh, X>::Type
    tanh(const ScalarClosure<X> &x);

//-- acos(x) -------------------------------------------------------------------
template <typename X>
    const typename _SUO_<OpArccos, X>::Type
    acos(const ScalarClosure<X> &x);

//-- asin(x) -------------------------------------------------------------------
template <typename X>
    const typename _SUO_<OpArcsin, X>::Type
    asin(const ScalarClosure<X> &x);

//-- atan(x) -------------------------------------------------------------------
template <typename X>
    const typename _SUO_<OpArctan, X>::Type
    atan(const ScalarClosure<X> &x);

} // namespace flens

#include <flens/scalaroperations.tcc>

#endif // FLENS_SCALAROPERATIONS_H
