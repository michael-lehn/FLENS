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

#include <cmath>
#include <complex>

namespace flens {

//== auxiliary traits for ScalarClosures =======================================

//-- ScalarOperation -----------------------------------------------------------

template <typename Op, typename X, typename Y>
struct _SO_
{
    typedef ScalarOperation<Op, X, Y>  Engine;
    typedef ScalarClosure<Engine>      Type;
};

//-- ScalarUnaryOperation ------------------------------------------------------

template <typename Op, typename X>
struct _SUO_
{
    typedef ScalarUnaryOperation<Op, X>  Engine;
    typedef ScalarClosure<Engine>        Type;
};

//== binary operators ==========================================================

//-- x + y ---------------------------------------------------------------------
template <typename X, typename Y>
const typename _SO_<OpAdd, X, Y>::Type
operator+(const ScalarClosure<X> &x, const ScalarClosure<Y> &y)
{
    typedef typename _SO_<OpAdd, X, Y>::Engine SO;
    return SO(x.engine(), y.engine());
}

template <typename X, typename T>
const typename _SO_<OpAdd, X, Scalar<T> >::Type
operator+(const ScalarClosure<X> &x, const T &a)
{
    typedef typename _SO_<OpAdd, X, Scalar<T> >::Engine SO;
    return SO(x.engine(), a);
}

template <typename T, typename X>
const typename _SO_<OpAdd, Scalar<T>, X>::Type
operator+(const T &a, const ScalarClosure<X> &x)
{
    typedef typename _SO_<OpAdd, Scalar<T>, X>::Engine SO;
    return SO(a, x.engine());
}

template <>
struct Operation<OpAdd>
{
    template <typename X, typename Y>
    static typename Promotion<X, Y>::Type
    eval(const X &x, const Y &y)
    {
        return x+y;
    }
};

//-- x - y ---------------------------------------------------------------------
// see below: implementation of unary operator -x

//-- x * y ---------------------------------------------------------------------
template <typename X, typename Y>
const typename _SO_<OpMult, X, Y>::Type
operator*(const ScalarClosure<X> &x, const ScalarClosure<Y> &y)
{
    typedef typename _SO_<OpMult, X, Y>::Engine SO;
    return SO(x.engine(), y.engine());
}

template <typename X, typename T>
const typename _SO_<OpMult, X, Scalar<T> >::Type
operator*(const ScalarClosure<X> &x, const T &a)
{
    typedef typename _SO_<OpMult, X, Scalar<T> >::Engine SO;
    return SO(x.engine(), a);
}

template <typename T, typename X>
const typename _SO_<OpMult, Scalar<T>, X>::Type
operator*(const T &a, const ScalarClosure<X> &x)
{
    typedef typename _SO_<OpMult, Scalar<T>, X>::Engine SO;
    return SO(a, x.engine());
}

template <>
struct Operation<OpMult>
{
    template <typename X, typename Y>
    static typename Promotion<X, Y>::Type
    eval(const X &x, const Y &y)
    {
        return x*y;
    }
};

//-- x / y----------------------------------------------------------------------
template <typename X, typename Y>
const typename _SO_<OpDiv, X, Y>::Type
operator/(const ScalarClosure<X> &x, const ScalarClosure<Y> &y)
{
    typedef typename _SO_<OpDiv, X, Y>::Engine SO;
    return SO(x.engine(), y.engine());
}

template <typename X, typename T>
const typename _SO_<OpDiv, X, Scalar<T> >::Type
operator/(const ScalarClosure<X> &x, const T &a)
{
    typedef typename _SO_<OpDiv, X, Scalar<T> >::Engine SO;
    return SO(x.engine(), a);
}

template <typename T, typename X>
const typename _SO_<OpDiv, Scalar<T>, X>::Type
operator/(const T &a, const ScalarClosure<X> &x)
{
    typedef typename _SO_<OpDiv, Scalar<T>, X>::Engine SO;
    return SO(a, x.engine());
}

template <>
struct Operation<OpDiv>
{
    template <typename X, typename Y>
    static typename Promotion<X, Y>::Type
    eval(const X &x, const Y &y)
    {
        return x/y;
    }
};

//-- x % y----------------------------------------------------------------------
template <typename X, typename Y>
const typename _SO_<OpMod, X, Y>::Type
operator%(const ScalarClosure<X> &x, const ScalarClosure<Y> &y)
{
    typedef typename _SO_<OpMod, X, Y>::Engine SO;
    return SO(x.engine(), y.engine());
}

template <typename X, typename T>
const typename _SO_<OpMod, X, Scalar<T> >::Type
operator%(const ScalarClosure<X> &x, const T &a)
{
    typedef typename _SO_<OpMod, X, Scalar<T> >::Engine SO;
    return SO(x.engine(), a);
}

template <typename T, typename X>
const typename _SO_<OpMod, Scalar<T>, X>::Type
operator%(const T &a, const ScalarClosure<X> &x)
{
    typedef typename _SO_<OpMod, Scalar<T>, X>::Engine SO;
    return SO(a, x.engine());
}

template <>
struct Operation<OpMod>
{
    template <typename X, typename Y>
    static typename Promotion<X, Y>::Type
    eval(const X &x, const Y &y)
    {
        return x%y;
    }
};

//-- max(x,y) ------------------------------------------------------------------
template <typename X, typename Y>
const typename _SO_<OpMax, X, Y>::Type
max(const ScalarClosure<X> &x, const ScalarClosure<Y> &y)
{
    typedef typename _SO_<OpMax, X, Y>::Engine SO;
    return SO(x.engine(), y.engine());
}

template <typename X>
const typename _SO_<OpMax, X, X>::Type
max(const ScalarClosure<X> &x, const ScalarClosure<X> &y)
{
    typedef typename _SO_<OpMax, X, X>::Engine SO;
    return SO(x.engine(), y.engine());
}

template <typename X, typename T>
const typename _SO_<OpMax, X, Scalar<T> >::Type
max(const ScalarClosure<X> &x, const T &a)
{
    typedef typename _SO_<OpMax, X, Scalar<T> >::Engine SO;
    return SO(x.engine(), a);
}

template <typename T, typename X>
const typename _SO_<OpMax, Scalar<T>, X>::Type
max(const T &a, const ScalarClosure<X> &x)
{
    typedef typename _SO_<OpMax, Scalar<T>, X>::Engine SO;
    return SO(a, x.engine());
}

template <>
struct Operation<OpMax>
{
    template <typename X>
    static const X &
    eval(const X &x, const X &y)
    {
        return std::max(x,y);
    }

    template <typename X, typename Y>
    static typename Promotion<X, Y>::Type
    eval(const X &x, const Y &y)
    {
        typedef typename Promotion<X, Y>::Type T;
        const T &X = x;
        const T &Y = y;
        return std::max(X,Y);
    }
};

//-- min(x,y) ------------------------------------------------------------------
template <typename X, typename Y>
const typename _SO_<OpMin, X, Y>::Type
min(const ScalarClosure<X> &x, const ScalarClosure<Y> &y)
{
    typedef typename _SO_<OpMin, X, Y>::Engine SO;
    return SO(x.engine(), y.engine());
}

template <typename X>
const typename _SO_<OpMin, X, X>::Type
min(const ScalarClosure<X> &x, const ScalarClosure<X> &y)
{
    typedef typename _SO_<OpMin, X, X>::Engine SO;
    return SO(x.engine(), y.engine());
}

template <typename X, typename T>
const typename _SO_<OpMin, X, Scalar<T> >::Type
min(const ScalarClosure<X> &x, const T &a)
{
    typedef typename _SO_<OpMin, X, Scalar<T> >::Engine SO;
    return SO(x.engine(), a);
}

template <typename T, typename X>
const typename _SO_<OpMin, Scalar<T>, X>::Type
min(const T &a, const ScalarClosure<X> &x)
{
    typedef typename _SO_<OpMin, Scalar<T>, X>::Engine SO;
    return SO(a, x.engine());
}

template <>
struct Operation<OpMin>
{
    template <typename X>
    static const X &
    eval(const X &x, const X &y)
    {
        return std::min(x,y);
    }

    template <typename X, typename Y>
    static typename Promotion<X, Y>::Type
    eval(const X &x, const Y &y)
    {
        typedef typename Promotion<X, Y>::Type T;
        const T &X = x;
        const T &Y = y;
        return std::min(X,Y);
    }
};

//-- pow(x,y) ------------------------------------------------------------------
template <typename X, typename Y>
const typename _SO_<OpPow, X, Y>::Type
pow(const ScalarClosure<X> &x, const ScalarClosure<Y> &y)
{
    typedef typename _SO_<OpPow, X, Y>::Engine SO;
    return SO(x.engine(), y.engine());
}

template <typename X>
const typename _SO_<OpPow, X, X>::Type
pow(const ScalarClosure<X> &x, const ScalarClosure<X> &y)
{
    typedef typename _SO_<OpPow, X, X>::Engine SO;
    return SO(x.engine(), y.engine());
}

template <typename X, typename T>
const typename _SO_<OpPow, X, Scalar<T> >::Type
pow(const ScalarClosure<X> &x, const T &a)
{
    typedef typename _SO_<OpPow, X, Scalar<T> >::Engine SO;
    return SO(x.engine(), a);
}

template <typename T, typename X>
const typename _SO_<OpPow, Scalar<T>, X>::Type
pow(const T &a, const ScalarClosure<X> &x)
{
    typedef typename _SO_<OpPow, Scalar<T>, X>::Engine SO;
    return SO(a, x.engine());
}

template <>
struct Operation<OpPow>
{
    template <typename X, typename Y>
    static typename Promotion<X, Y>::Type
    eval(const X &x, const Y &y)
    {
        return std::pow(x,y);
    }
};

//== unary and binary operators ================================================

//-- x - y, -x -----------------------------------------------------------------
template <typename X, typename Y>
const typename _SO_<OpMinus, X, Y>::Type
operator-(const ScalarClosure<X> &x, const ScalarClosure<Y> &y)
{
    typedef typename _SO_<OpMinus, X, Y>::Engine SO;
    return SO(x.engine(), y.engine());
}

template <typename X, typename T>
const typename _SO_<OpMinus, X, Scalar<T> >::Type
operator-(const ScalarClosure<X> &x, const T &a)
{
    typedef typename _SO_<OpMinus, X, Scalar<T> >::Engine SO;
    return SO(x.engine(), a);
}

template <typename T, typename X>
const typename _SO_<OpMinus, Scalar<T>, X>::Type
operator-(const T &a, const ScalarClosure<X> &x)
{
    typedef typename _SO_<OpMinus, Scalar<T>, X>::Engine SO;
    return SO(a, x.engine());
}

template <typename X>
const typename _SUO_<OpMinus, X>::Type
operator-(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpMinus, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpMinus>
{
    // binary operation: x - y
    template <typename X, typename Y>
    static typename Promotion<X, Y>::Type
    eval(const X &x, const Y &y)
    {
        return x-y;
    }

    // unary operation: -x
    template <typename X>
    static X
    eval(const X &x)
    {
        return -x;
    }
};

//== unary operators ==========================================================

//-- conjugate(x) or conj(x)----------------------------------------------------
template <typename X>
const typename _SUO_<OpConj, X>::Type
conjugate(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpConj, X>::Engine SUO;
    return SUO(x.engine());
}

template <typename X>
const typename _SUO_<OpConj, X>::Type
conj(const ScalarClosure<X> &x)
{
    const typename _SUO_<OpConj, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpConj>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::conj(x);
    }
};

//-- ceil ----------------------------------------------------------------------
template <typename X>
const typename _SUO_<OpCeil, X>::Type
ceil(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpCeil, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpCeil>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::ceil(x);
    }
};

//-- floor ---------------------------------------------------------------------
template <typename X>
const typename _SUO_<OpFloor, X>::Type
floor(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpFloor, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpFloor>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::floor(x);
    }
};

//-- abs(x) --------------------------------------------------------------------
template <typename X>
const typename _SUO_<OpAbs, X>::Type
abs(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpAbs, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpAbs>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::abs(x);
    }
};

//-- sqrt(x) -------------------------------------------------------------------
template <typename X>
const typename _SUO_<OpSqrt, X>::Type
sqrt(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpSqrt, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpSqrt>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::sqrt(x);
    }
};

//-- exp(x) --------------------------------------------------------------------
template <typename X>
const typename _SUO_<OpExp, X>::Type
exp(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpExp, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpExp>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::exp(x);
    }
};

//-- log(x) --------------------------------------------------------------------
template <typename X>
const typename _SUO_<OpLog, X>::Type
log(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpLog, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpLog>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::log(x);
    }
};

//-- cos(x) --------------------------------------------------------------------
template <typename X>
const typename _SUO_<OpCos, X>::Type
cos(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpCos, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpCos>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::cos(x);
    }
};

//-- sin(x) --------------------------------------------------------------------
template <typename X>
const typename _SUO_<OpSin, X>::Type
sin(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpSin, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpSin>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::sin(x);
    }
};

//-- tan(x) --------------------------------------------------------------------
template <typename X>
const typename _SUO_<OpTan, X>::Type
tan(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpTan, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpTan>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::tan(x);
    }
};

//-- cosh(x) -------------------------------------------------------------------
template <typename X>
const typename _SUO_<OpCosh, X>::Type
cosh(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpCosh, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpCosh>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::cosh(x);
    }
};

//-- sinh(x) -------------------------------------------------------------------
template <typename X>
const typename _SUO_<OpSinh, X>::Type
sinh(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpSinh, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpSinh>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::sinh(x);
    }
};

//-- tanh(x) -------------------------------------------------------------------
template <typename X>
const typename _SUO_<OpTanh, X>::Type
tanh(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpTanh, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpTanh>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::tanh(x);
    }
};

//-- acos(x) -------------------------------------------------------------------
template <typename X>
const typename _SUO_<OpArccos, X>::Type
acos(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpArccos, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpArccos>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::acos(x);
    }
};

//-- asin(x) -------------------------------------------------------------------
template <typename X>
const typename _SUO_<OpArcsin, X>::Type
asin(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpArcsin, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpArcsin>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::asin(x);
    }
};

//-- atan(x) -------------------------------------------------------------------
template <typename X>
const typename _SUO_<OpArctan, X>::Type
atan(const ScalarClosure<X> &x)
{
    typedef typename _SUO_<OpArctan, X>::Engine SUO;
    return SUO(x.engine());
}

template <>
struct Operation<OpArctan>
{
    template <typename X>
    static X
    eval(const X &x)
    {
        return std::atan(x);
    }
};

} // namespace flens
