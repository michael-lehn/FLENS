/*
 *   Copyright (c) 2012, Michael Lehn, Klaus Pototzky
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

#ifndef PLAYGROUND_FLENS_MPI_BCAST_BCAST_TCC
#define PLAYGROUND_FLENS_MPI_BCAST_BCAST_TCC 1

#include<playground/flens/mpi/mpi-flens.h>

namespace flens { namespace mpi {

template <typename VX>
typename RestrictTo<IsDenseVector<VX>::value,
                    void>::Type
MPI_bcast(VX &&x, const int root = 0)
{
#ifdef WITH_MPI
    using namespace MPI;
    
    typedef typename RemoveRef<VX>::Type   VectorX;
    typedef typename VectorX::ElementType T;
    typedef typename VectorX::IndexType   IndexType;
    
    IndexType length = x.length();
    
    COMM_WORLD.Bcast(&length, MPI_Type<IndexType>::size, MPI_Type<IndexType>::Type, root); 
    
    if ( x.length()== 0) {
        x.resize(length);
    }
    
    ASSERT( x.length()==length );
    
    for (IndexType i=x.firstIndex(); i<=x.lastIndex(); ++i) {
             COMM_WORLD.Bcast(reinterpret_cast<typename MPI_Type<T>::PrimitiveType *>(&x(i)), 
                              MPI_Type<T>::size, MPI_Type<T>::Type, root);  
    }
 
#endif
}
  
  
template <typename MA>
typename RestrictTo<IsGeMatrix<MA>::value,
                    void>::Type
MPI_bcast(MA &&A, const int root = 0)
{
#ifdef WITH_MPI
    using namespace MPI;
    
    typedef typename RemoveRef<MA>::Type   MatrixA;
    typedef typename MatrixA::IndexType    IndexType;
    
    const Underscore<IndexType> _;

    IndexType numRows = A.numRows();
    IndexType numCols = A.numCols();
    
    COMM_WORLD.Bcast(&numRows, MPI_Type<IndexType>::size, MPI_Type<IndexType>::Type, root); 
    COMM_WORLD.Bcast(&numCols, MPI_Type<IndexType>::size, MPI_Type<IndexType>::Type, root); 
    
    if ( A.numRows()==0 && A.numCols()==0 ) {
        A.resize(numRows, numCols); 
    }
    
    ASSERT( A.numRows()==numRows );
    ASSERT( A.numCols()==numCols );
    
    for (IndexType i=A.firstCol(); i<=A.lastCol(); ++i) {
        auto x = A(_,i); 
        MPI_bcast(x, root);
    }
#endif
}

} }

#endif // PLAYGROUND_FLENS_MPI_BCAST_BCAST_TCC
