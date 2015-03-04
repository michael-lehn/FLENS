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

#ifndef FLENS_BLAS_LEVEL1EXTENSIONS_ASUM1_TCC
#define FLENS_BLAS_LEVEL1EXTENSIONS_ASUM1_TCC 1

#include <cxxblas/cxxblas.h>
#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/closures/closures.h>
#include <flens/blas/level1/level1.h>
#include <flens/typedefs.h>


namespace flens { namespace blas {

//-- BLAS Level 1 --------------------------------------------------------------

template <typename X, typename T>
typename RestrictTo<IsNotComplex<T>::value, void>::Type
asum1(const DenseVector<X> &x, T &absoluteSum)
{
#   ifdef HAVE_CXXBLAS_ASUM1
    cxxblas::asum1(x.length(), x.data(), x.stride(), absoluteSum);
#   else
    ASSERT(0);
#   endif
}

template <typename X>
const typename ComplexTrait<typename X::ElementType>::PrimitiveType
asum1(const DenseVector<X> &x)
{
    typename ComplexTrait<typename X::ElementType>::PrimitiveType  absoluteSum;

    asum1(x, absoluteSum);
    return absoluteSum;
}

//-- BLAS Level 1 extensions ---------------------------------------------------

//== GeneralMatrix

template <typename MA, typename T>
typename RestrictTo<IsNotComplex<T>::value, void>::Type
asum1(const GeMatrix<MA> &A, T &absoluteSum)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    const Underscore<IndexType>  _;

    absoluteSum = T();
    if (A.order()==ColMajor) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            absoluteSum += asum1(A(_,j));
        }
    } else {
        for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
            absoluteSum += asum1(A(i,_));
        }
    }
}

template <typename MA>
const typename ComplexTrait<typename MA::ElementType>::PrimitiveType
asum1(const GeMatrix<MA> &A)
{
    typename ComplexTrait<typename MA::ElementType>::PrimitiveType  absoluteSum;

    asum1(A, absoluteSum);
    return absoluteSum;
}


} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1EXTENSIONS_ASUM1_TCC
