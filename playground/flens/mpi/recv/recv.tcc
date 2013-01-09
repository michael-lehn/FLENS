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

#ifdef WITH_MPI

template <typename T>
typename RestrictTo<MPI_Type<T>::Compatible,
                        T>::Type
MPI_recv(const int source, const MPI::Comm &communicator)
{
    T val(0);
     
    MPI_recv(val, source, communicator);

    return val;
}


template <typename T>
typename RestrictTo<MPI_Type<T>::Compatible,
                    void>::Type
MPI_recv(T &x, const int source, const MPI::Comm &communicator)
{

    using namespace MPI;
    
    Status status;
    communicator.Recv(reinterpret_cast<typename MPI_Type<T>::PrimitiveType *>(&x), 
                      MPI_Type<T>::size, MPI_Type<T>::Type, source, 0, status);
        
    return ;


}

template <typename IndexType, typename T>
typename RestrictTo<MPI_Type<T>::Compatible,
                    void>::Type
MPI_recv(const IndexType n, T *x, const IndexType incX, const int source,
         const MPI::Comm &communicator)
{

    using namespace MPI;

    
    if ( incX==IndexType(1) ) {
        Status status;
        communicator.Recv(reinterpret_cast<typename MPI_Type<T>::PrimitiveType *>(x), 
                          n*MPI_Type<T>::size, MPI_Type<T>::Type, source, 0, status);
    } else {
        
        for (IndexType i=0, iX=0; i<n; ++i, iX+=incX) {
             MPI_recv(x[iX], source, communicator);
        }
    }
    return ;

}


template <typename VX>
typename RestrictTo<IsDenseVector<VX>::value,
                    void>::Type
MPI_recv(VX &&x, const int source, const MPI::Comm &communicator)
{

    using namespace MPI;
    
    typedef typename RemoveRef<VX>::Type   VectorX;
    typedef typename VectorX::ElementType T;
    typedef typename VectorX::IndexType   IndexType;
    

    IndexType length, stride;
        
    // Receive Vector length and stride
    MPI_recv(length, source, communicator);
    MPI_recv(stride, source, communicator);
    
    if ( x.length()==0 ) {
        x.resize(length); 
    }
    
    ASSERT (x.length()==length );
    ASSERT (x.stride()==stride || x.stride() != 1);
    
    MPI_recv(x.length(), x.data(), x.stride(), source, communicator);

}

template <typename MA>
typename RestrictTo<IsGeMatrix<MA>::value,
                    void>::Type
MPI_recv(MA &&A, const int source, const MPI::Comm &communicator)
{

    using namespace MPI;
    
    typedef typename RemoveRef<MA>::Type   MatrixA;
    typedef typename MatrixA::ElementType T;
    typedef typename MatrixA::IndexType   IndexType;
    
    
    IndexType numCols, numRows;
        
    // Receive size
    MPI_recv(numRows, source, communicator);
    MPI_recv(numCols, source, communicator);  
    
    if ( A.numCols()==0 && A.numRows()==0 ) {
        A.resize(numRows, numCols);
    }
    
    ASSERT ( A.numRows()==numRows );
    ASSERT ( A.numCols()==numCols );
    
#ifndef NDEBUG
    int isColMajor;
    MPI_recv(isColMajor, source, communicator);
    
    ASSERT( isColMajor == ( A.order() == ColMajor ) );
#endif    
    
    if ( A.order() == ColMajor ) {
      
      for (IndexType i=0; i<A.numCols(); ++i) {

            MPI_recv(A.numRows(), A.data()+i*A.leadingDimension(), 1, source, communicator);
        
        }
      
    } else {
      
        for (IndexType i=0; i<A.numRows(); ++i) {

            MPI_recv(A.numCols(), A.data()+i*A.leadingDimension(), 1, source, communicator);
        
        
        }      
    }


}

template <typename X>
typename RestrictTo<IsDenseVector<X>::value || 
                    IsGeMatrix<X>::value,
                    void>::Type
MPI_recv_all(std::vector<X> &x, const MPI::Comm &communicator)
{

    if (x.size()==0) {
        x.resize(MPI_size());   
    }
    ASSERT( x.size() == MPI_size() );
    
    for (int i=0; i<MPI_size(); ++i) {
        MPI_recv(x.at(i), i, communicator);
    }
    
}

#endif // WITH_MPI

} }

#endif // PLAYGROUND_FLENS_MPI_RECV_RECV_TCC
