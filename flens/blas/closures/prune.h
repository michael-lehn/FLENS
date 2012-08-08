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

#ifndef FLENS_BLAS_CLOSURES_PRUNE_H
#define FLENS_BLAS_CLOSURES_PRUNE_H 1

#include <cxxblas/typedefs.h>
#include <flens/auxiliary/promotion.h>
#include <flens/blas/operators/operators.h>
#include <flens/matrixtypes/impl/matrixclosure.h>
#include <flens/typedefs.h>

namespace flens {

//== PruneScaling ==============================================================

//-- General definition

template <typename MatVec>
struct PruneScaling
{
    typedef MatVec                                  Remainder;
    typedef typename Remainder::Impl::ElementType   ScalingType;

    static const ScalingType
    getFactor(const MatVec &)
    {
        return ScalingType(1);
    }

    static const Remainder &
    getRemainder(const MatVec &x)
    {
        return x;
    }
};

//-- Specialization for scaled matrix or vector

template <typename T, typename M>
struct PruneScaling<MatrixClosure<OpMult, ScalarValue<T>, M> >
{
    typedef MatrixClosure<OpMult, ScalarValue<T>, M>  MC;
    typedef M                                         Remainder;
    typedef T                                         ScalingType;

    static const ScalingType
    getFactor(const MC &x)
    {
        return x.left().value();
    }

    static const Remainder &
    getRemainder(const MC &x)
    {
        return x.right();
    }
};

template <typename T, typename V>
struct PruneScaling<VectorClosure<OpMult, ScalarValue<T>, V> >
{
    typedef VectorClosure<OpMult, ScalarValue<T>, V>  VC;
    typedef V                                         Remainder;
    typedef T                                         ScalingType;

    static const ScalingType
    getFactor(const VC &x)
    {
        return x.left().value();
    }

    static const Remainder &
    getRemainder(const VC &x)
    {
        return x.right();
    }
};

//== PruneConjTrans ============================================================

//-- General definition

template <typename Matrix>
struct PruneConjTrans
{
    typedef Matrix Remainder;

    static const Transpose trans = NoTrans;

    static const Remainder &
    remainder(const Matrix &matrix)
    {
        return matrix;
    }

};

//-- Specialization for OpTrans, OpConj closures

template <typename M>
struct PruneConjTrans<MatrixClosure<OpTrans, M, M> >
{
    typedef MatrixClosure<OpTrans, M, M>           MC;
    typedef typename PruneConjTrans<M>::Remainder  Remainder;

    static const Transpose trans = Transpose(Trans^PruneConjTrans<M>::trans);

    static const Remainder &
    remainder(const MC &A)
    {
        return PruneConjTrans<M>::remainder(A.left());
    }

};

template <typename M>
struct PruneConjTrans<MatrixClosure<OpConj, M, M> >
{
    typedef MatrixClosure<OpConj, M, M>            MC;
    typedef typename PruneConjTrans<M>::Remainder  Remainder;

    static const Transpose trans = Transpose(Conj^PruneConjTrans<M>::trans);

    static const Remainder &
    remainder(const MC &A)
    {
        return PruneConjTrans<M>::remainder(A.left());
    }
};

} // namespace flens

#endif // FLENS_BLAS_CLOSURES_PRUNE_H
