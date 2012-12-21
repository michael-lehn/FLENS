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

#ifndef PLAYGROUND_FLENS_MPI_SEND_SEND_TCC
#define PLAYGROUND_FLENS_MPI_SEND_SEND_TCC 1

#include<playground/flens/mpi/mpi-flens.h>

namespace flens { namespace mpi {

template <typename VX>
typename RestrictTo<IsDenseVector<VX>::value,
                    void>::Type
MPI_send(VX &&x, const int source, const int dest)
{
#ifdef WITH_MPI
    using namespace MPI;
    
    typedef typename RemoveRef<VX>::Type   VectorX;
    typedef typename VectorX::ElementType T;
    typedef typename VectorX::IndexType   IndexType;
    

    if ( source != MPI_rank()) {
        return; 
    }
    
    // Send Vector length and stride
    const IndexType length = x.length();
    const IndexType stride = x.stride();
    COMM_WORLD.Send(&length, MPI_Type<IndexType>::size, MPI_Type<IndexType>::Type, dest, 0);
    COMM_WORLD.Send(&stride, MPI_Type<IndexType>::size, MPI_Type<IndexType>::Type, dest, 0);    
        
    if (x.stride()==1) {
      COMM_WORLD.Send(reinterpret_cast<typename MPI_Type<T>::PrimitiveType *>(x.data()), 
                          x.length()*MPI_Type<T>::size, 
                          MPI_Type<T>::Type, dest, 0);
    } else {
        // Send Data
        for (IndexType i=x.firstIndex(); i<=x.lastIndex(); ++i) {
          COMM_WORLD.Send(reinterpret_cast<typename MPI_Type<T>::PrimitiveType *>(&x(i)), 
                          MPI_Type<T>::size, 
                          MPI_Type<T>::Type, dest, 0);
            
        }
    }
#else
    ASSERT(0);
#endif
}

template <typename MA>
typename RestrictTo<IsGeMatrix<MA>::value,
                    void>::Type
MPI_send(MA &&A, const int source, const int dest)
{
#ifdef WITH_MPI
    using namespace MPI;
    
    typedef typename RemoveRef<MA>::Type   MatrixA;
    typedef typename MatrixA::ElementType T;
    typedef typename MatrixA::IndexType   IndexType;
    

    if ( source != MPI_rank()) {
        return; 
    }
    
    // Send size
    const IndexType numCols = A.numCols();
    const IndexType numRows = A.numRows();
    
    COMM_WORLD.Send(&numRows, MPI_Type<IndexType>::size, MPI_Type<IndexType>::Type, dest, 0);
    COMM_WORLD.Send(&numCols, MPI_Type<IndexType>::size, MPI_Type<IndexType>::Type, dest, 0);  

#ifndef NDEBUG
    const int isColMajor = ( A.order() == ColMajor );
    COMM_WORLD.Send(&isColMajor, MPI_Type<int>::size, MPI_Type<int>::Type, dest, 0);
#endif
    
    if ( A.order() == ColMajor ) {
      
        for (IndexType i=0; i<A.numCols(); ++i) {

            COMM_WORLD.Send(reinterpret_cast<typename MPI_Type<T>::PrimitiveType *>(A.data()+i*A.leadingDimension()), 
                      A.numRows()*MPI_Type<T>::size, 
                      MPI_Type<T>::Type, dest, 0);
        
        
        }
      
    } else {
      
        for (IndexType i=0; i<A.numRows(); ++i) {

            COMM_WORLD.Send(reinterpret_cast<typename MPI_Type<T>::PrimitiveType *>(A.data()+i*A.leadingDimension()), 
                      A.numCols()*MPI_Type<T>::size, 
                      MPI_Type<T>::Type, dest, 0);
        
        
        }      
    }

#else
    ASSERT(0);
#endif
}


template <typename VX>
typename RestrictTo<IsDenseVector<VX>::value,
                    void>::Type
MPI_send(VX &&x, const int dest)
{
#ifdef WITH_MPI  
    for (int i=0; i<MPI_size(); ++i) {
        MPI_send(x, i, dest); 
    }
#else
    ASSERT(0);    
#endif
}

template <typename MA>
typename RestrictTo<IsGeMatrix<MA>::value,
                    void>::Type
MPI_send(MA &&A, const int dest)
{
#ifdef WITH_MPI  
    for (int i=0; i<MPI_size(); ++i) {
        MPI_send(A, i, dest); 
    }
#else
    ASSERT(0);    
#endif
}

} }

#endif // PLAYGROUND_FLENS_MPI_SEND_SEND_TCC
