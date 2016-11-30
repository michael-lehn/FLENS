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

#ifndef FLENS_BLAS_OPERATORS_OPCROSS_H
#define FLENS_BLAS_OPERATORS_OPCROSS_H 1

#include <flens/vectortypes/vectortypes.h>

namespace flens {

struct OpCross {};

//-- cross product -------------------------------------------------------------
// (x % y)
template <typename VX, typename VY>
    const VectorClosure<OpCross, typename VX::Impl, typename VY::Impl>
    operator%(const Vector<VX> &x, const Vector<VY> &y);

//-- triple vector products ----------------------------------------------------
// x * (y % z)
template <typename VX, typename VY, typename VZ>
    typename Promotion<typename VX::Impl::ElementType,
                       typename Promotion<typename VY::Impl::ElementType,
                                          typename VZ::Impl::ElementType>::Type>::Type
    operator*(const VX &x, const VectorClosure<OpCross, VY, VZ> &yz);

// (x % y) * z
template <typename VX, typename VY, typename VZ>
    typename Promotion<typename Promotion<typename VX::Impl::ElementType,
                                          typename VY::Impl::ElementType>::Type,
                       typename VZ::Impl::ElementType>::Type
    operator*(const VectorClosure<OpCross, VX, VY> &xy, const VZ &z);

} // namespace flens

#endif // FLENS_BLAS_OPERATORS_OPCROSS_H
