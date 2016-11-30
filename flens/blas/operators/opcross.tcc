/*
 *   Copyright (c) 2016, Thibaud Kloczko
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

#ifndef FLENS_BLAS_OPERATORS_OPCROSS_TCC
#define FLENS_BLAS_OPERATORS_OPCROSS_TCC 1

#include <flens/blas/operators/opcross.h>

namespace flens {

//-- cross product -------------------------------------------------------------
// (x % y)

template <typename VX, typename VY>
const VectorClosure<OpCross, typename VX::Impl, typename VY::Impl>
operator%(const Vector<VX> &x, const Vector<VY> &y)
{
    typedef VectorClosure<OpCross, typename VX::Impl, typename VY::Impl>  VC;
    return VC(x.impl(), y.impl());
}

//-- triple vector products ----------------------------------------------------
// x * (y % z)
template <typename VX, typename VY, typename VZ>
    typename Promotion<typename VX::Impl::ElementType,
                       typename Promotion<typename VY::Impl::ElementType,
                                          typename VZ::Impl::ElementType>::Type>::Type
    operator*(const VX &x, const VectorClosure<OpCross, VY, VZ> &yz)
{
    typedef typename VX::Impl::ElementType TX;
    typedef typename VY::Impl::ElementType TY;
    typedef typename VZ::Impl::ElementType TZ;

    typedef typename Promotion<TX, typename Promotion<TY, TZ>::Type>::Type T;

    const typename Result<VX>::Type& _x = x;
    const typename Result<VY>::Type& _y = yz.left();
    const typename Result<VZ>::Type& _z = yz.right();

    T res = _x(1) * (_y(2) * _z(3) - _y(3) * _z(2));
    res  += _x(2) * (_y(3) * _z(1) - _y(1) * _z(3));
    res  += _x(3) * (_y(1) * _z(2) - _y(2) * _z(1));

    return res;
}

// (x % y) * z
template <typename VX, typename VY, typename VZ>
    typename Promotion<typename Promotion<typename VX::Impl::ElementType,
                                          typename VY::Impl::ElementType>::Type,
                       typename VZ::Impl::ElementType>::Type
    operator*(const VectorClosure<OpCross, VX, VY> &xy, const VZ &z)
{
    typedef VectorClosure<OpCross, VX, VY> VC;

    const typename Result<VC>::Type xy_ = xy;
    const typename Result<VZ>::Type& z_ = z;

    return (xy_ * z_);
}

} // namespace flens

#endif // FLENS_BLAS_OPERATORS_OPCROSS_TCC
