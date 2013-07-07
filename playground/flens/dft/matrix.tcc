/*
 *   Copyright (c) 2013, Klaus Pototzky
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

#ifndef PLAYGROUND_FLENS_DFT_MATRIX_TCC
#define PLAYGROUND_FLENS_DFT_MATRIX_TCC 1

#include<playground/cxxdft/direction.h>

namespace flens {
namespace dft {

template <typename AIN, typename AOUT>
typename RestrictTo<IsComplexGeMatrix<AIN>::value &&
                    IsComplexGeMatrix<AOUT>::value,
                    void>::Type
dft_col_forward(AIN &&Ain, AOUT &&Aout)
{

    typedef typename RemoveRef<AOUT>::Type          MatrixA;
    typedef typename MatrixA::IndexType             IndexType;
    
    
    ASSERT( Ain.numCols()==Aout.numCols() );
    ASSERT( Ain.numRows()==Aout.numRows() );

    const IndexType numRows = Ain.numRows();
    const IndexType numCols = Ain.numCols();

    const IndexType AinStride = ( Ain.order()==StorageOrder::ColMajor ? 1 : Ain.leadingDimension() );
    const IndexType AinDist   = ( Ain.order()==StorageOrder::RowMajor ? 1 : Ain.leadingDimension() );

    const IndexType AoutStride = ( Aout.order()==StorageOrder::ColMajor ? 1 : Aout.leadingDimension() );
    const IndexType AoutDist   = ( Aout.order()==StorageOrder::RowMajor ? 1 : Aout.leadingDimension() );
    
    cxxdft::dft_multiple(numRows, numCols, 
                         Ain.data(), AinStride, AinDist,
                         Aout.data(), AoutStride, AoutDist,
                         cxxdft::DFTDirection::Forward);

}

template <typename AIN, typename AOUT>
typename RestrictTo<IsComplexGeMatrix<AIN>::value &&
                    IsComplexGeMatrix<AOUT>::value,
                    void>::Type
dft_col_backward(AIN &&Ain, AOUT &&Aout)
{

    typedef typename RemoveRef<AOUT>::Type          MatrixA;
    typedef typename MatrixA::IndexType             IndexType;
    
    ASSERT( Ain.numCols()==Aout.numCols() );
    ASSERT( Ain.numRows()==Aout.numRows() );

    const IndexType numRows = Ain.numRows();
    const IndexType numCols = Ain.numCols();

    const IndexType AinStride = ( Ain.order()==StorageOrder::ColMajor ? 1 : Ain.leadingDimension() );
    const IndexType AinDist   = ( Ain.order()==StorageOrder::RowMajor ? 1 : Ain.leadingDimension() );

    const IndexType AoutStride = ( Aout.order()==StorageOrder::ColMajor ? 1 : Aout.leadingDimension() );
    const IndexType AoutDist   = ( Aout.order()==StorageOrder::RowMajor ? 1 : Aout.leadingDimension() );
    
    cxxdft::dft_multiple(numRows, numCols, 
                         Ain.data(), AinStride, AinDist,
                         Aout.data(), AoutStride, AoutDist,
                         cxxdft::DFTDirection::Backward);

}
    
template <typename AIN, typename AOUT>
typename RestrictTo<IsComplexGeMatrix<AIN>::value &&
                    IsComplexGeMatrix<AOUT>::value,
                    void>::Type
dft_row_forward(AIN &&Ain, AOUT &&Aout)
{

    typedef typename RemoveRef<AOUT>::Type          MatrixA;
    typedef typename MatrixA::IndexType             IndexType;
    
    ASSERT( Ain.numCols()==Aout.numCols() );
    ASSERT( Ain.numRows()==Aout.numRows() );

    const IndexType numRows = Ain.numRows();
    const IndexType numCols = Ain.numCols();

    const IndexType AinStride = ( Ain.order()==StorageOrder::RowMajor ? 1 : Ain.leadingDimension() );
    const IndexType AinDist   = ( Ain.order()==StorageOrder::ColMajor ? 1 : Ain.leadingDimension() );

    const IndexType AoutStride = ( Aout.order()==StorageOrder::RowMajor ? 1 : Aout.leadingDimension() );
    const IndexType AoutDist   = ( Aout.order()==StorageOrder::ColMajor ? 1 : Aout.leadingDimension() );
    
    cxxdft::dft_multiple(numCols, numRows, 
                         Ain.data(), AinStride, AinDist,
                         Aout.data(), AoutStride, AoutDist,
                         cxxdft::DFTDirection::Forward);

}

template <typename AIN, typename AOUT>
typename RestrictTo<IsComplexGeMatrix<AIN>::value &&
                    IsComplexGeMatrix<AOUT>::value,
                    void>::Type
dft_row_backward(AIN &&Ain, AOUT &&Aout)
{

    typedef typename RemoveRef<AOUT>::Type          MatrixA;
    typedef typename MatrixA::IndexType             IndexType;
    
    ASSERT( Ain.numCols()==Aout.numCols() );
    ASSERT( Ain.numRows()==Aout.numRows() );

    const IndexType numRows = Ain.numRows();
    const IndexType numCols = Ain.numCols();

    const IndexType AinStride = ( Ain.order()==StorageOrder::RowMajor ? 1 : Ain.leadingDimension() );
    const IndexType AinDist   = ( Ain.order()==StorageOrder::ColMajor ? 1 : Ain.leadingDimension() );

    const IndexType AoutStride = ( Aout.order()==StorageOrder::RowMajor ? 1 : Aout.leadingDimension() );
    const IndexType AoutDist   = ( Aout.order()==StorageOrder::ColMajor ? 1 : Aout.leadingDimension() );
    
    cxxdft::dft_multiple(numCols, numRows, 
                         Ain.data(), AinStride, AinDist,
                         Aout.data(), AoutStride, AoutDist,
                         cxxdft::DFTDirection::Backward);

}



template <typename AIN, typename AOUT>
typename RestrictTo<IsComplexGeMatrix<AIN>::value &&
                    IsComplexGeMatrix<AOUT>::value,
                    void>::Type
dft_col_forward_normalized(AIN &&Ain, AOUT &&Aout)
{
    
    typedef typename RemoveRef<AOUT>::Type          MatrixA;
    typedef typename MatrixA::ElementType           T;
    typedef typename ComplexTrait<T>::PrimitiveType PT;
    
    dft_col_forward(Ain, Aout);
    Aout /= PT(Ain.numRows());
    
}

template <typename AIN, typename AOUT>
typename RestrictTo<IsComplexGeMatrix<AIN>::value &&
                    IsComplexGeMatrix<AOUT>::value,
                    void>::Type
dft_col_backward_normalized(AIN &&Ain, AOUT &&Aout)
{
    
    typedef typename RemoveRef<AOUT>::Type          MatrixA;
    typedef typename MatrixA::ElementType           T;
    typedef typename ComplexTrait<T>::PrimitiveType PT;
    
    dft_col_backward(Ain, Aout);
    Aout /= PT(Ain.numRows());
    
}

template <typename AIN, typename AOUT>
typename RestrictTo<IsComplexGeMatrix<AIN>::value &&
                    IsComplexGeMatrix<AOUT>::value,
                    void>::Type
dft_row_forward_normalized(AIN &&Ain, AOUT &&Aout)
{
    
    typedef typename RemoveRef<AOUT>::Type          MatrixA;
    typedef typename MatrixA::ElementType           T;
    typedef typename ComplexTrait<T>::PrimitiveType PT;
    
    dft_row_forward(Ain, Aout);
    Aout /= PT(Ain.numCols());
    
}

template <typename AIN, typename AOUT>
typename RestrictTo<IsComplexGeMatrix<AIN>::value &&
                    IsComplexGeMatrix<AOUT>::value,
                    void>::Type
dft_row_backward_normalized(AIN &&Ain, AOUT &&Aout)
{
    
    typedef typename RemoveRef<AOUT>::Type          MatrixA;
    typedef typename MatrixA::ElementType           T;
    typedef typename ComplexTrait<T>::PrimitiveType PT;
    
    dft_row_backward(Ain, Aout);
    Aout /= PT(Ain.numCols());
    
}
    
} // namespace dft
} // namespace flens

#endif // PLAYGROUND_FLENS_DFT_MATRIX_TCC
