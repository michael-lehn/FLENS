/*
 *   Copyright (c) 2012, Michael Lehn
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

#ifndef FLENS_BLAS_CLOSURES_CLOSURETYPE_H
#define FLENS_BLAS_CLOSURES_CLOSURETYPE_H 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/closures/prune.h>
#include <flens/blas/operators/operators.h>

namespace flens { namespace blas {

template <typename Op, typename L, typename R>
struct ClosureType
{
    static const bool isMatrixVectorProduct = false;
    static const bool isResidual            = false;
    static const bool isMatrixMatrixProduct = false;
};

//
//  A*x (matrix-vector product)
//  x*A (vector-matrix product)
//  A*B (matrix-matrix product)
//
template <typename L, typename R>
struct ClosureType<OpMult, L, R>
{
    static const bool isMatrixVectorProduct = IsMatrix<L>::value
                                           && IsVector<R>::value;

    static const bool isVectorMatrixProduct = IsVector<L>::value
                                           && IsMatrix<R>::value;

    static const bool isMatrixMatrixProduct = IsMatrix<L>::value
                                           && IsMatrix<R>::value;
};

//
//  y + A*x  or  beta*y + A*x  (matrix-vector product)
//
template <typename VY, typename MA, typename VX>
struct ClosureType<OpAdd, VY, VectorClosure<OpMult, MA, VX> >
{
    // remove scaling
    typedef typename PruneScaling<VY>::Remainder   RVY;

    static const bool vectorY  = IsVector<VY>::value
                             && !IsClosure<RVY>::value;

    static const bool isMatrixVectorProduct = vectorY
                                           && IsMatrix<MA>::value
                                           && IsVector<VX>::value;
};

//
//  C + A*B  or  beta*C + A*B  (matrix-matrix product)
//
template <typename MC, typename MA, typename MB>
struct ClosureType<OpAdd, MC, MatrixClosure<OpMult, MA, MB> >
{
    // remove scaling
    typedef typename PruneScaling<MC>::Remainder   RMC;

    static const bool matrixC  = IsMatrix<MC>::value
                             && !IsClosure<RMC>::value;

    static const bool isMatrixMatrixProduct = matrixC
                                           && IsMatrix<MA>::value
                                           && IsMatrix<MB>::value;
};

//
//  b - A*x
//
template <typename VB, typename MA, typename VX>
struct ClosureType<OpSub, VB, VectorClosure<OpMult, MA, VX> >
{
    static const bool isMatrixVectorProduct = IsVector<VB>::value
                                           && IsMatrix<MA>::value
                                           && IsVector<VX>::value;

    static const bool isResidual = isMatrixVectorProduct
                                && !IsClosure<VB>::value
                                && !IsClosure<MA>::value
                                && !IsClosure<VX>::value;
};


} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_CLOSURETYPE_H
