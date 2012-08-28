/*
 *   Copyright (c) 2010, Michael Lehn
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

#ifndef FLENS_BLAS_CLOSURES_PRUNEMATRIXCLOSURE_H
#define FLENS_BLAS_CLOSURES_PRUNEMATRIXCLOSURE_H 1

#include <cxxblas/typedefs.h>
#include <flens/auxiliary/constref.h>
#include <flens/auxiliary/promotion.h>
#include <flens/blas/operators/operators.h>
#include <flens/matrixtypes/impl/matrixclosure.h>
#include <flens/typedefs.h>

namespace flens {

//-- General definition --------------------------------------------------------
template <typename Matrix>
struct PruneMatrixClosure
{
    typedef int     ScalingFactor;
    typedef Matrix  Remainder;

    template <typename ALPHA>
        static const ALPHA &
        updateScalingFactor(const ALPHA &alpha, const Matrix &);

    static Transpose
    updateTranspose(Transpose trans);

    static typename ConstRef<Remainder>::Type
    remainder(const Matrix &matrix);
};

//-- Specialization for particular closures ------------------------------------
template <typename L, typename R>
    struct PruneMatrixClosure<MatrixClosure<OpMult, ScalarValue<L>, R> >;

template <typename R>
    struct PruneMatrixClosure<MatrixClosure<OpTrans, R, R> >;

template <typename R>
    struct PruneMatrixClosure<MatrixClosure<OpConj, R, R> >;

} // namespace flens

#endif // FLENS_BLAS_CLOSURES_PRUNEMATRIXCLOSURE_H
