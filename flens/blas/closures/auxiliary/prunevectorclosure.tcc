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

#ifndef FLENS_BLAS_CLOSURES_PRUNEVECTORCLOSURE_TCC
#define FLENS_BLAS_CLOSURES_PRUNEVECTORCLOSURE_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/closures/auxiliary/prunevectorclosure.h>
#include <flens/typedefs.h>

namespace flens {

//-- General definition --------------------------------------------------------

template <typename Vector>
template <typename ALPHA>
const ALPHA &
PruneVectorClosure<Vector>::updateScalingFactor(const ALPHA &alpha,
                                                const Vector &)
{
    return alpha;
}

template <typename Vector>
Transpose
PruneVectorClosure<Vector>::updateTranspose(Transpose trans)
{
    return trans;
}

template <typename Vector>
typename ConstRef<typename PruneVectorClosure<Vector>::Remainder>::Type
PruneVectorClosure<Vector>::remainder(const Vector &vector)
{
    return vector;
}

//-- Specialization for particular closures ------------------------------------
//-- Closure from alpha*A
template <typename L, typename R>
struct PruneVectorClosure<VectorClosure<OpMult, ScalarValue<L>, R> >
{
    typedef VectorClosure<OpMult, ScalarValue<L>, R>  VC;

    typedef typename PruneVectorClosure<R>::ScalingFactor  _ScalingFactor;
    typedef typename Promotion<L, _ScalingFactor>::Type    ScalingFactor;
    typedef typename PruneVectorClosure<R>::Remainder      Remainder;

    template <typename ALPHA>
    static const typename Promotion<ALPHA, ScalingFactor>::Type
    updateScalingFactor(const ALPHA &alpha, const VC &vc)
    {
        typedef PruneVectorClosure<R> PVC;
        return PVC::updateScalingFactor(alpha*vc.left().value(), vc.right());
    }

    static Transpose
    updateTranspose(Transpose trans)
    {
        return PruneVectorClosure<R>::updateTranspose(trans);
    }

    static typename ConstRef<Remainder>::Type
    remainder(const VC &vc)
    {
        return PruneVectorClosure<R>::remainder(vc.right());
    }
};

//-- Closure from conj(A)
template <typename R>
struct PruneVectorClosure<VectorClosure<OpConj, R, R> >
{
    typedef VectorClosure<OpConj, R, R>  VC;

    typedef typename PruneVectorClosure<R>::ScalingFactor  ScalingFactor;
    typedef typename PruneVectorClosure<R>::Remainder      Remainder;

    template <typename ALPHA>
    static const typename Promotion<ALPHA, ScalingFactor>::Type
    updateScalingFactor(const ALPHA &alpha, const VC &vc)
    {
        return PruneVectorClosure<R>::updateScalingFactor(alpha, vc.right());
    }

    static Transpose
    updateTranspose(Transpose trans)
    {
        trans = PruneVectorClosure<R>::updateTranspose(trans);
        return Transpose(trans^Conj);
    }

    static typename ConstRef<Remainder>::Type
    remainder(const VC &vc)
    {
        return PruneVectorClosure<R>::remainder(vc.right());
    }
};

} // namespace flens

#endif // FLENS_BLAS_CLOSURES_PRUNEVECTORCLOSURE_TCC
