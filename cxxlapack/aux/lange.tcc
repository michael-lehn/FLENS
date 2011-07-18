/*
 *   Copyright (c) 2011, Iris Haecker
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


/*
 *   Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.
 *   
 *   $COPYRIGHT$
 *   
 *   Additional copyrights may follow
 *   
 *   $HEADER$
 *   
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions are
 *   met:
 *   
 *   - Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer. 
 *     
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer listed
 *     in this license in the documentation and/or other materials
 *     provided with the distribution.
 *     
 *   - Neither the name of the copyright holders nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
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

#ifndef CXXLAPACK_AUX_LANGE_TCC
#define CXXLAPACK_AUX_LANGE_TCC 1

#include <cxxblas/cxxblas.h>
#include <cxxlapack/cxxlapack.h>

namespace cxxlapack {

using cxxblas::asum;

template <typename IndexType, typename MA, typename T>
void
lange(StorageOrder order, Norm norm,
      IndexType m, IndexType n,
      const MA *A, IndexType ldA,
      T *work, T &value)
{
    assert(order==ColMajor);

    if (min(m,n)==0) {
        value = 0;
    } else if (norm==MaximumNorm) {
//
//      Find max(abs(A(i,j))).
//
        value = 0;
        for (IndexType j=0; j<n; ++j) {
            for (IndexType i=0; i<m; ++i) {
                value = max(value, abs(A[i+ldA*j]));
            }
        }
    } else if (norm==OneNorm) {
//
//      Find norm1(A).
//
        value = 0;
        for (IndexType j=0; j<n; ++j) {
            T sum;

            asum(m, A+j*ldA, IndexType(1), sum);
            value = max(value, sum);
        }
    } else if (norm==InfinityNorm) {
//
//      Find normI(A).
//
        for (IndexType i=0; i<m; ++i) {
            work[i] = 0;
        }
        for (IndexType j=0; j<n; ++j) {
            for (IndexType i=0; i<m; ++i) {
                work[i] += abs(A[i+ldA*j]);
            }
        }
        value = 0;
        for (IndexType i=0; i<m; ++i) {
            value = max(value, work[i]);
        }
    } else if (norm==FrobeniusNorm) {
//
//      Find normF(A).
//
        T scale = 0;
        T sum = 1;
        for (IndexType j=0; j<n; ++j) {
            lassq(m, A+j*ldA, IndexType(1), scale, sum);
        }
        value = scale*sqrt(sum);
    }
}

} // namespace cxxlapack

#endif // CXXLAPACK_AUX_LANGE_TCC
