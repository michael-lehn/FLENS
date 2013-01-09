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

#ifdef WITH_MPI  
template <typename T>
typename RestrictTo<MPI_Type<T>::Compatible,
                    void>::Type
MPI_send(const T &x, const int dest, const MPI::Comm &communicator)
{

    using namespace MPI;
    
    communicator.Send(
             reinterpret_cast<const typename MPI_Type<T>::PrimitiveType *>(&x),
             MPI_Type<T>::size, MPI_Type<T>::Type, dest, 0);
        
    return ;
    
}

template <typename IndexType, typename T>
typename RestrictTo<MPI_Type<T>::Compatible,
                    void>::Type
MPI_send(const IndexType n, const T *x, const IndexType incX, const int dest, 
         const MPI::Comm &communicator)
{

    using namespace MPI;
    
    
    if ( incX==1 ) {
        communicator.Send(
                 reinterpret_cast<const typename MPI_Type<T>::PrimitiveType *>(x),
                 n*MPI_Type<T>::size, MPI_Type<T>::Type, dest, 0);
    } else {
        for (IndexType i=0, iX=0; i<n; ++i, iX+=incX) {
	    MPI_send(x[iX], dest, communicator); 
	}
    }
    
    return ;
    
}

template <typename VX>
typename RestrictTo<IsDenseVector<VX>::value,
                    void>::Type
MPI_send(VX &&x, const int dest, const MPI::Comm &communicator)
{
   using namespace MPI;
    
    typedef typename RemoveRef<VX>::Type   VectorX;
    typedef typename VectorX::ElementType T;
    typedef typename VectorX::IndexType   IndexType;
    
    
    // Send Vector length and stride
    const IndexType length = x.length();
    const IndexType stride = x.stride();   
    
    MPI_send(length, dest, communicator);
    MPI_send(stride, dest, communicator); 
    
    MPI_send(x.length(), x.data(), x.stride(), dest, communicator);
    
}

template <typename MA>
typename RestrictTo<IsGeMatrix<MA>::value,
                    void>::Type
MPI_send(MA &&A, const int dest, const MPI::Comm &communicator)
{

    using namespace MPI;
    
    typedef typename RemoveRef<MA>::Type   MatrixA;
    typedef typename MatrixA::ElementType T;
    typedef typename MatrixA::IndexType   IndexType;
    
    // Send size
    const IndexType numCols = A.numCols();
    const IndexType numRows = A.numRows();
    
    MPI_send(numRows, dest, communicator);
    MPI_send(numCols, dest, communicator);

#ifndef NDEBUG
    const int isColMajor = ( A.order() == ColMajor );
    MPI_send(isColMajor, dest, communicator);
#endif
    
    if ( A.order() == ColMajor ) {
      
        for (IndexType i=0; i<A.numCols(); ++i) {
            MPI_send(A.numRows(), A.data()+i*A.leadingDimension(), 1,
                     dest, communicator);     
        }
      
    } else {
      
        for (IndexType i=0; i<A.numRows(); ++i) {

            MPI_send(A.numCols(), A.data()+i*A.leadingDimension(), 1,
                     dest, communicator);
        
        
        }      
    }

}

#endif // WITH_MPI

} }

#endif // PLAYGROUND_FLENS_MPI_SEND_SEND_TCC
