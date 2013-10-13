/*
 *   Copyright (c) 2012, Klaus Pototzky
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
 *
 */

#ifndef PLAYGROUND_FLENS_LAPACKEXTENSIONS_SY_DETERMINANT_TCC
#define PLAYGROUND_FLENS_LAPACKEXTENSIONS_SY_DETERMINANT_TCC 1

#include <playground/flens/blas-extensions/blas-extensions.h>
#include <playground/flens/lapack-extensions/sy/determinant.h>

namespace flens { namespace lapack { namespace extensions {

//-- det(sy)
template <typename MA, typename VPIV, typename VWORK>
typename RestrictTo< IsSyMatrix<MA>::value
                  && IsIntegerDenseVector<VPIV>::value
                  && IsDenseVector<VWORK>::value,
typename RemoveRef<MA>::Type::ElementType>::Type
det(MA &&A, VPIV &&piv, VWORK &&work)
{
    ASSERT(A.numCols()==A.numRows());

    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::ElementType   T;
    typedef typename MatrixA::IndexType     IndexType;

    trf(A, piv, work);
    T value(1);

    for (IndexType i=A.firstRow(), k=piv.firstIndex(); i<=A.lastRow(); ++i, ++k) {
        if (piv(k)>0) {
            value *= A(i,i);
        } else if (i < A.lastRow() && piv(k)<0 && piv(k)==piv(k+1) ) {
            if (A.upLo()==Upper) {
                value *= (A(i+1,i+1)*A(i,i)-A(i,i+1)*A(i,i+1));
            } else {
                value *= (A(i+1,i+1)*A(i,i)-A(i+1,i)*A(i+1,i));
            }
            ++k;
            ++i;
        }

    }
    return value;
}

template <typename MA, typename VPIV>
typename RestrictTo< IsSyMatrix<MA>::value
                  && IsIntegerDenseVector<VPIV>::value,
typename RemoveRef<MA>::Type::ElementType>::Type
det(MA &&A, VPIV &&piv)
{
    typedef typename RemoveRef<MA>::Type::Vector Vector;

    Vector  work;
    return det(A, piv, work);
}

template <typename MA>
typename RestrictTo<IsSyMatrix<MA>::value,
typename RemoveRef<MA>::Type::ElementType>::Type
det(MA &&A)
{
    typedef typename RemoveRef<MA>::Type         MatrixA;
    typedef typename MatrixA::IndexType          IndexType;
    typedef typename MatrixA::Vector             Vector;

    Vector                          work;
    DenseVector<Array<IndexType> >  piv;

    return det(A, piv, work);
}

} } } // namespace extensions, lapack, flens

#endif // PLAYGROUND_FLENS_LAPACKEXTENSIONS_SY_DETERMINANT_TCC
