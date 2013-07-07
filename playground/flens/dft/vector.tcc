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

#ifndef PLAYGROUND_FLENS_DFT_VECTOR_TCC
#define PLAYGROUND_FLENS_DFT_VECTOR_TCC 1

#include<playground/cxxdft/direction.h>

namespace flens {
namespace dft {
    
template <typename VIN, typename VOUT>
typename RestrictTo<IsComplexDenseVector<VIN>::value &&
                    IsComplexDenseVector<VOUT>::value,
                    void>::Type
dft_forward(VIN &&vin, VOUT &&vout)
{
	//
	//  Remove references from rvalue types
	//
	typedef typename RemoveRef<VIN>::Type    VectorV;
	typedef typename VectorV::IndexType      IndexType;
    
    ASSERT( vin.length()==vout.length() );
    
    const IndexType N = vin.length();
    
    cxxdft::dft_single(N, 
                       vin.data(), vin.stride(), 
                       vout.data(), vout.stride(), 
                       cxxdft::DFTDirection::Forward);

}

template <typename VIN, typename VOUT>
typename RestrictTo<IsComplexDenseVector<VIN>::value &&
                    IsComplexDenseVector<VOUT>::value,
                    void>::Type
dft_backward(VIN &&vin, VOUT &&vout)
{
	//
	//  Remove references from rvalue types
	//
	typedef typename RemoveRef<VIN>::Type    VectorV;
	typedef typename VectorV::IndexType      IndexType;
    
    ASSERT( vin.length()==vout.length() );
    
    const IndexType N = vin.length();
    
    cxxdft::dft_single(N, 
                       vin.data(), vin.stride(), 
                       vout.data(), vout.stride(), 
                       cxxdft::DFTDirection::Backward);

}

template <typename VIN, typename VOUT>
typename RestrictTo<IsComplexDenseVector<VIN>::value &&
                    IsComplexDenseVector<VOUT>::value,
                    void>::Type
dft_forward_normalized(VIN &&vin, VOUT &&vout)
{
    
    typedef typename RemoveRef<VOUT>::Type          VectorV;
    typedef typename VectorV::ElementType           T;
    typedef typename ComplexTrait<T>::PrimitiveType PT;
    
    dft_forward(vin, vout);
    vout /= PT(vin.length());
    
}

template <typename VIN, typename VOUT>
typename RestrictTo<IsComplexDenseVector<VIN>::value &&
                    IsComplexDenseVector<VOUT>::value,
                    void>::Type
dft_backward_normalized(VIN &&vin, VOUT &&vout)
{
    
    typedef typename RemoveRef<VOUT>::Type          VectorV;
    typedef typename VectorV::ElementType           T;
    typedef typename ComplexTrait<T>::PrimitiveType PT;
    
    dft_backward(vin, vout);
    vout /= PT(vin.length());
    
}
    
} // namespace dft
} // namespace flens

#endif // PLAYGROUND_FLENS_DFT_VECTOR_TCC
