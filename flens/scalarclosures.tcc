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

#include <flens/scalaroperations.h>

namespace flens {

// == ScalarClosure ============================================================


template <typename Engine>
ScalarClosure<Engine>::ScalarClosure()
{
}

template <typename Engine>
ScalarClosure<Engine>::ScalarClosure(const Engine &engine)
    : _engine(engine)
{
}

template <typename Engine>
ScalarClosure<Engine>::~ScalarClosure()
{
}

template <typename Engine>
void
ScalarClosure<Engine>::operator=(const ScalarClosure<Engine> &rhs)
{
    copy(rhs, engine());
}

template <typename Engine>
template <typename RHS>
void
ScalarClosure<Engine>::operator=(const ScalarClosure<RHS> &rhs)
{
    copy(rhs, engine());
}

template <typename Engine>
const Engine &
ScalarClosure<Engine>::engine() const
{
    return _engine;
}

template <typename Engine>
Engine &
ScalarClosure<Engine>::engine()
{
    return _engine;
}

//== Scalar ====================================================================

template <typename T>
Scalar<T>::Scalar(T value)
    : _value(value)
{
}

template <typename T>
Scalar<T>::operator T() const
{
    return _value;
}

template <typename T>
const T&
Scalar<T>::value() const
{
    return _value;
}

// == ScalarOperation ==========================================================

template <typename Op, typename L, typename R>
ScalarOperation<Op,L,R>::ScalarOperation(const L &l, const R &r)
    : _left(l), _right(r)
{
}

template <typename Op, typename L, typename R>
typename ScalarOperation<Op,L,R>::T
ScalarOperation<Op,L,R>::value() const
{
    return Operation<Op>::eval(_left.value(), _right.value());
}

// == ScalarUnaryOperation =====================================================

template <typename Op, typename Arg>
ScalarUnaryOperation<Op, Arg>::ScalarUnaryOperation(const Arg &arg)
    : _arg(arg)
{
}

template <typename Op, typename Arg>
typename ScalarUnaryOperation<Op, Arg>::T
ScalarUnaryOperation<Op, Arg>::value() const
{
    return Operation<Op>::eval(_arg.value());
}

// == VectorElement ============================================================

template <typename VecImpl>
VectorElement<VecImpl>::VectorElement(VecImpl &vecImpl, Index &index)
    : _vecImpl(vecImpl), _index(index)
{
}

template <typename VecImpl>
const typename VectorElement<VecImpl>::T &
VectorElement<VecImpl>::value() const
{
    return _vecImpl(_index.engine().value());
}

//------------------------------------------------------------------------------

template <typename X, typename Y>
void
copy(const ScalarClosure<X> &x, VectorElement<Y> &_y)
{
    typedef typename Y::ElementType T;

    Y             &y = _y._vecImpl;
    Index::Engine &i = _y._index.engine();

    for (i.value()=y.firstIndex(); i.value()<=y.lastIndex(); ++i.value()) {
        y(i.value()) = x.engine().value();
    }
}

// == MatrixElement ============================================================

template <typename MatImpl>
MatrixElement<MatImpl>::MatrixElement(MatImpl &matImpl, Index &i, Index &j)
    : _matImpl(matImpl), _i(i), _j(j)
{
}

template <typename MatImpl>
const typename MatrixElement<MatImpl>::T &
MatrixElement<MatImpl>::value() const
{
    return _matImpl(_i.engine().value(), _j.engine().value());
}

//------------------------------------------------------------------------------

template <typename X, typename Y>
void
copy(const ScalarClosure<X> &x, MatrixElement<Y> &_y)
{
    typedef typename Y::ElementType T;

    Y             &y = _y._matImpl;
    Index::Engine &i = _y._i.engine();
    Index::Engine &j = _y._j.engine();

    // loop (kind of) optimized for colmajor storage
    // TODO: specialize for rowmajor
    for (j.value()=y.firstCol(); j.value()<=y.lastCol(); ++(j.value())) {
        for (i.value()=y.firstRow(); i.value()<=y.lastRow(); ++(i.value())) {
            y(i.value(),j.value()) = x.engine().value();
        }
    }
}

// == Summation ================================================================

template <typename I>
Summation<I>::Summation(const ScalarClosure<I> &summand, Index &index)
    : _summand(summand.engine()), _index(index)
{
}

template <typename I>
typename Summation<I>::T
Summation<I>::value() const
{
    Index::Engine &i = _index.engine();
    T             value = 0;
    for (i.value()=i.firstIndex(); i.value()<=i.lastIndex(); ++i.value()) {
        value += _summand.value();
    }
    return value;
}

//------------------------------------------------------------------------------

template <typename I>
const ScalarClosure<Summation<I> >
sum(const ScalarClosure<I> &summand, Index &index)
{
    return Summation<I>(summand, index);
}

//==============================================================================

} // namespace flens
