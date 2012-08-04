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

#ifndef FLENS_BLAS_OPERATORS_OPMINUS_TCC
#define FLENS_BLAS_OPERATORS_OPMINUS_TCC 1

#include <flens/blas/closures/closures.h>
#include <flens/blas/level1/dot.h>
#include <flens/blas/operators/opminus.h>
#include <flens/typedefs.h>

namespace flens {

template <typename VX>
const VectorClosure<OpMult,
                     ScalarValue<typename VX::Impl::ElementType>,
                     typename VX::Impl>
operator-(const Vector<VX> &x)
{
    typedef typename VX::Impl::ElementType    T;
    typedef VectorClosure<OpMult, ScalarValue<T>, typename VX::Impl>  VC;

    return VC(T(-1), x.impl());
}

template <typename MA>
const MatrixClosure<OpMult,
                    ScalarValue<typename MA::Impl::ElementType>,
                    typename MA::Impl>
operator-(const Matrix<MA> &A)
{
    typedef typename MA::Impl::ElementType    T;
    typedef MatrixClosure<OpMult, ScalarValue<T>, typename MA::Impl>  MC;

    return MC(T(-1), A.impl());
}

} // namespace flens

#endif // FLENS_BLAS_OPERATORS_OPMINUS_TCC
