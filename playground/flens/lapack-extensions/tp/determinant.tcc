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

#ifndef PLAYGROUND_FLENS_LAPACKEXTENSIONS_TP_DETERMINANT_TCC
#define PLAYGROUND_FLENS_LAPACKEXTENSIONS_TP_DETERMINANT_TCC 1

namespace flens { namespace lapack { namespace extensions {

//-- det(tb)
template <typename MA>
typename RestrictTo<IsTpMatrix<MA>::value,
typename RemoveRef<MA>::Type::ElementType>::Type
det(MA &&A)
{
    
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::ElementType   T;
    typedef typename MatrixA::IndexType     IndexType;
    
    if (A.diag()==Unit)
        return T(1);
    
    T value(1);
    
    for (IndexType i=A.firstRow(), j=A.firstCol();i<=A.lastRow();++i, ++j) {
        value *= A(i,j);
    }

    return value;
}

} } } // namespace extensions, lapack, flens

#endif // PLAYGROUND_FLENS_LAPACKEXTENSIONS_TP_DETERMINANT_TCC
