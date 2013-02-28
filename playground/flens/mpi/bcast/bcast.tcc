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

#ifdef WITH_MPI

template <typename T>
typename RestrictTo<MPI_Type<T>::Compatible,
                    void>::Type
MPI_bcast(T &x, const int root, const MPI::Comm &communicator)
{
    communicator.Bcast(reinterpret_cast<typename MPI_Type<T>::PrimitiveType *>(&x), 
                     MPI_Type<T>::size, MPI_Type<T>::Type(), root);
}
  
template <typename IndexType, typename T>
typename RestrictTo<MPI_Type<T>::Compatible,
                    void>::Type
MPI_bcast(const IndexType n, T *x, const IndexType incX, const int root, 
          const MPI::Comm &communicator)
{
  
    if ( incX==1 ) {
        communicator.Bcast(reinterpret_cast<typename MPI_Type<T>::PrimitiveType *>(x), 
                           n*MPI_Type<T>::size, MPI_Type<T>::Type(), root);
    } else {
    
        for (IndexType i=0, iX=0; i<n; ++i, iX+=incX) {
            MPI_bcast(x[iX], root, communicator);
        }
    
    }

}

template <typename VX>
typename RestrictTo<IsDenseVector<VX>::value,
                    void>::Type
MPI_bcast(VX &&x, const int root, const MPI::Comm &communicator)
{


    
    typedef typename RemoveRef<VX>::Type   VectorX;
    typedef typename VectorX::ElementType T;
    typedef typename VectorX::IndexType   IndexType;
    
    IndexType length = x.length();
    IndexType stride = x.stride();    
    MPI_bcast(length, root, communicator); 
    MPI_bcast(stride, root, communicator); 
    
    if ( x.length()== 0) {
        x.resize(length);
    }
    
    ASSERT( x.length()==length );
    
    ASSERT( x.stride()==stride || stride!=1 );
    MPI_bcast(x.length(), x.data(), x.stride(), root, communicator);
 

}
  
  
template <typename MA>
typename RestrictTo<IsGeMatrix<MA>::value,
                    void>::Type
MPI_bcast(MA &&A, const int root, const MPI::Comm &communicator)
{

    
    typedef typename RemoveRef<MA>::Type   MatrixA;
    typedef typename MatrixA::IndexType    IndexType;
    
    const Underscore<IndexType> _;

    IndexType numRows = A.numRows();
    IndexType numCols = A.numCols();
    
    MPI_bcast(numRows, root, communicator); 
    MPI_bcast(numCols, root, communicator); 
    
    if ( A.numRows()==0 && A.numCols()==0 ) {
        A.resize(numRows, numCols); 
    }
    
    ASSERT( A.numRows()==numRows );
    ASSERT( A.numCols()==numCols );
    
    if ( A.order() == ColMajor ) {
        for (IndexType i=A.firstCol(); i<=A.lastCol(); ++i) {
            auto x = A(_,i); 
            MPI_bcast(x.length(), x.data(), x.stride(), root, communicator);
        }
    } else {
        for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
            auto x = A(i,_); 
            MPI_bcast(x.length(), x.data(), x.stride(), root, communicator);
        }
    }

}


#else

template <typename T>
typename RestrictTo<(IsInteger<T>::value ||
                     IsReal<T>::value ||
                     IsComplex<T>::value),
                    void>::Type
MPI_bcast(T &x, const int root)
{
    ASSERT( root==0 );  
}

template <typename T>
typename RestrictTo<(IsDenseVector<T>::value ||
                     IsGeMatrix<T>::value) ,
                     void>::Type
MPI_bcast(T &&x, const int root)
{
    ASSERT( root==0 );  
}

#endif

} }

#endif // PLAYGROUND_FLENS_MPI_BCAST_BCAST_TCC
