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

#ifndef PLAYGROUND_FLENS_MPI_RECV_RECV_TCC
#define PLAYGROUND_FLENS_MPI_RECV_RECV_TCC 1

#include<playground/flens/mpi/mpi-flens.h>

namespace flens { namespace mpi {

template <typename VX>
typename RestrictTo<IsDenseVector<VX>::value,
                    void>::Type
MPI_recv(VX &&x, const int source, const int dest)
{
#ifdef WITH_MPI
    using namespace MPI;
    
    typedef typename RemoveRef<VX>::Type   VectorX;
    typedef typename VectorX::ElementType T;
    typedef typename VectorX::IndexType   IndexType;
    
    if ( dest != MPI_rank()) {
        return; 
    }
    
    IndexType length, stride;
    
    Status status;
    
    // Receive Vector length and stride
    
    COMM_WORLD.Recv(&length, MPI_Type<IndexType>::size, MPI_Type<IndexType>::Type, source, 0, status);
    COMM_WORLD.Recv(&stride, MPI_Type<IndexType>::size, MPI_Type<IndexType>::Type, source, 0, status);
    
    if ( x.length()==0 ) {
        x.resize(length); 
    }
    
    ASSERT (x.length()==length );
    
    if ( stride==1) {
      
      ASSERT( x.stride()==stride );
      
      COMM_WORLD.Recv(reinterpret_cast<typename MPI_Type<T>::PrimitiveType *>(x.data()), 
                            x.length()*MPI_Type<T>::size, MPI_Type<T>::Type, source, 0, status);
    } else {
        // Receive Data 
        for (IndexType i=x.firstIndex(); i<=x.lastIndex(); ++i) {
            COMM_WORLD.Recv(reinterpret_cast<typename MPI_Type<T>::PrimitiveType *>(&x(i)), 
                            MPI_Type<T>::size, MPI_Type<T>::Type, source, 0, status);
        }
    }
#else
    ASSERT(0);    
#endif
}

template <typename MA>
typename RestrictTo<IsGeMatrix<MA>::value,
                    void>::Type
MPI_recv(MA &&A, const int source, const int dest)
{
#ifdef WITH_MPI
    using namespace MPI;
    
    typedef typename RemoveRef<MA>::Type   MatrixA;
    typedef typename MatrixA::ElementType T;
    typedef typename MatrixA::IndexType   IndexType;
    
    if ( dest != MPI_rank()) {
        return; 
    }
    
    IndexType numCols, numRows;
    
    Status status;
    
    // Receive size
    COMM_WORLD.Recv(&numRows, MPI_Type<IndexType>::size, MPI_Type<IndexType>::Type, source, 0, status);
    COMM_WORLD.Recv(&numCols, MPI_Type<IndexType>::size, MPI_Type<IndexType>::Type, source, 0, status);    
    
    if ( A.numCols()==0 && A.numRows()==0 ) {
        A.resize(numRows, numCols);
    }
    
    ASSERT ( A.numRows()==numRows );
    ASSERT ( A.numCols()==numCols );
    
#ifndef NDEBUG
    int isColMajor;
    COMM_WORLD.Recv(&isColMajor, MPI_Type<int>::size, MPI_Type<int>::Type, source, 0, status);
    
    ASSERT( isColMajor == ( A.order() == ColMajor ) );
#endif    
    
    if ( A.order() == ColMajor ) {
      
      for (IndexType i=0; i<A.numCols(); ++i) {

            COMM_WORLD.Recv(reinterpret_cast<typename MPI_Type<T>::PrimitiveType *>(A.data()+i*A.leadingDimension()), 
                      A.numRows()*MPI_Type<T>::size, 
                      MPI_Type<T>::Type, source, 0, status);
        
        
        }
      
    } else {
      
        for (IndexType i=0; i<A.numRows(); ++i) {

            COMM_WORLD.Recv(reinterpret_cast<typename MPI_Type<T>::PrimitiveType *>(A.data()+i*A.leadingDimension()), 
                      A.numCols()*MPI_Type<T>::size, 
                      MPI_Type<T>::Type, source, 0, status);
        
        
        }      
    }

#else
    ASSERT(0);
#endif
}

template <typename X>
typename RestrictTo<IsDenseVector<X>::value || 
                    IsGeMatrix<X>::value,
                    void>::Type
MPI_recv_all(std::vector<X> &x, const int root)
{
#ifdef WITH_MPI  
    if ( root != MPI_rank() ) {
        return;
    }
    if (x.size()==0) {
        x.resize(MPI_size());   
    }
    ASSERT( x.size() == MPI_size() );
    
    for (int i=0; i<MPI_size(); ++i) {
        MPI_recv(x.at(i), i, root);
    }
#else
    ASSERT(0);
#endif    
}



} }

#endif // PLAYGROUND_FLENS_MPI_RECV_RECV_TCC
