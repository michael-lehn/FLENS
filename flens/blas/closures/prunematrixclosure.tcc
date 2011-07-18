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

#ifndef FLENS_BLAS_CLOSURES_PRUNEMATRIXCLOSURE_TCC
#define FLENS_BLAS_CLOSURES_PRUNEMATRIXCLOSURE_TCC 1

namespace flens {

//-- General definition --------------------------------------------------------

template <typename Matrix>
template <typename ALPHA>
const ALPHA &
PruneMatrixClosure<Matrix>::updateScalingFactor(const ALPHA &alpha,
                                                const Matrix &)
{
    return alpha;
}

template <typename Matrix>
cxxblas::Transpose
PruneMatrixClosure<Matrix>::updateTranspose(cxxblas::Transpose trans)
{
    return trans;
}

template <typename Matrix>
const typename PruneMatrixClosure<Matrix>::Remainder &
PruneMatrixClosure<Matrix>::remainder(const Matrix &matrix)
{
    return matrix;
}

//-- Specialization for particular closures ------------------------------------ 
//-- Closure from alpha*A
template <typename L, typename R>
struct PruneMatrixClosure<MatrixClosure<OpMult, ScalarValue<L>, R> >
{
    typedef MatrixClosure<OpMult, ScalarValue<L>, R>  MC;
    
    typedef typename PruneMatrixClosure<R>::ScalingFactor  _ScalingFactor;
    typedef typename Promotion<L, _ScalingFactor>::Type    ScalingFactor;
    typedef typename PruneMatrixClosure<R>::Remainder      Remainder;
    
    template <typename ALPHA>
    static const typename Promotion<ALPHA, ScalingFactor>::Type
    updateScalingFactor(const ALPHA &alpha, const MC &mc)
    {
        typedef PruneMatrixClosure<R> PMC;
        return PMC::updateScalingFactor(alpha*mc.left().value(), mc.right());
    }

    static cxxblas::Transpose
    updateTranspose(cxxblas::Transpose trans)
    {
        return PruneMatrixClosure<R>::updateTranspose(trans);
    }

    static const Remainder &
    remainder(const MC &mc)
    {
        return PruneMatrixClosure<R>::remainder(mc.right());
    }
};

//-- Closure from transpose(A)
template <typename R>
struct PruneMatrixClosure<MatrixClosure<OpTrans, R, R> >
{
    typedef MatrixClosure<OpTrans, R, R>  MC;
    
    typedef typename PruneMatrixClosure<R>::ScalingFactor  ScalingFactor;
    typedef typename PruneMatrixClosure<R>::Remainder      Remainder;
    
    template <typename ALPHA>
    static const typename Promotion<ALPHA, ScalingFactor>::Type
    updateScalingFactor(const ALPHA &alpha, const MC &mc)
    {
        return PruneMatrixClosure<R>::updateScalingFactor(alpha, mc.right());
    }

    static cxxblas::Transpose
    updateTranspose(cxxblas::Transpose trans)
    {
        trans = PruneMatrixClosure<R>::updateTranspose(trans);
        return cxxblas::Transpose(trans^cxxblas::Trans);
    }

    static const Remainder &
    remainder(const MC &mc)
    {
        return PruneMatrixClosure<R>::remainder(mc.right());
    }
};

//-- Closure from conjugate(A)
template <typename R>
struct PruneMatrixClosure<MatrixClosure<OpConj, R, R> >
{
    typedef MatrixClosure<OpConj, R, R>  MC;
    
    typedef typename PruneMatrixClosure<R>::ScalingFactor  ScalingFactor;
    typedef typename PruneMatrixClosure<R>::Remainder      Remainder;
    
    template <typename ALPHA>
    static const typename Promotion<ALPHA, ScalingFactor>::Type
    updateScalingFactor(const ALPHA &alpha, const MC &mc)
    {
        return PruneMatrixClosure<R>::updateScalingFactor(alpha, mc.right());
    }

    static cxxblas::Transpose
    updateTranspose(cxxblas::Transpose trans)
    {
        trans = PruneMatrixClosure<R>::updateTranspose(trans);
        return cxxblas::Transpose(trans^cxxblas::Conj);
    }

    static const Remainder &
    remainder(const MC &mc)
    {
        return PruneMatrixClosure<R>::remainder(mc.right());
    }
};

} // namespace flens

#endif // FLENS_BLAS_CLOSURES_PRUNEMATRIXCLOSURE_TCC

