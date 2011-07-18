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

#ifndef CXXLAPACK_AUX_LAQGE_TCC
#define CXXLAPACK_AUX_LAQGE_TCC 1

#include <cxxblas/cxxblas.h>
#include <cxxlapack/cxxlapack.h>

namespace cxxlapack {

template <typename IndexType, typename T>
void
laqge(StorageOrder order,
      IndexType m, IndexType n, T *A, IndexType ldA,
      const T *r, const T *c, const T &rowCond, const T &colCond,
      const T &AMax, Equilibration &equilibration)
{
    assert(order==ColMajor);

    const T threshold = T(0.1);
//
//  Quick return if possible
//
    if ((m==0) || (n==0)) {
        equilibration = NoEquilibration;
        return;
    }
//
//  Initialize LARGE and SMALL.
//
    const T small = lamch<T>(SafeMin) / lamch<T>(Precision);
    const T large = T(1)/small;

    if ((rowCond>=threshold) && (AMax>=small) && (AMax<=large)) {
//
//      No row scaling
//
        if (colCond>=threshold) {
//
//          No column scaling
//
            equilibration = NoEquilibration;
        } else {
//
//          Column scaling
//
            for (IndexType j=0; j<n; ++j) {
                T cj = c[j];
                for (IndexType i=0; i<m; ++i) {
                    A[i+j*ldA] *= cj;
                }
            }
            equilibration = ColEquilibration;
        }
    } else if (colCond>=threshold) {
//
//      Row scaling, no column scaling
//
        for (IndexType j=0; j<n; ++j) {
            for (IndexType i=0; i<m; ++i) {
                A[i+j*ldA] *= r[i];
            }
        }
        equilibration = RowEquilibration;
    } else {
//
//      Row and column scaling
//
        for (IndexType j=0; j<n; ++j) {
            T cj = c[j];
            for (IndexType i=0; i<m; ++i) {
                A[i+j*ldA] *= cj*r[i];
            }
        }
        equilibration = RowColEquilibration;
    }
}

} // namespace cxxlapack

#endif // CXXLAPACK_AUX_LAQGE_TCC
