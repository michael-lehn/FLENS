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

#ifndef CXXLAPACK_GESV_GETF2_TCC
#define CXXLAPACK_GESV_GETF2_TCC 1

#include <algorithm>
#include <cxxblas/cxxblas.h>
#include <cxxlapack/cxxlapack.h>

namespace cxxlapack {

using cxxblas::ger;
using cxxblas::iamax;
using cxxblas::scal;
using cxxblas::swap;

template <typename IndexType, typename MA>
IndexType
getf2(StorageOrder order,
      IndexType m, IndexType n, MA *A, IndexType ldA,
      IndexType *iPiv)
{
    IndexType info = 0;

    if (m<0) {
        info = -1;
        // xerbla("getf2", -info);
        return info;
    }
    if (n<0) {
        info = -2;
        // xerbla("getf2", -info);
        return info;
    }
    if ((order==ColMajor)  && (ldA<max(IndexType(1),m))) {
        info = -4;
        // xerbla("getf2", -info);
        return info;
    }
//
//  Quick return if possible
//
    if ((m==0) || (n==0)) {
        return info;
    }
#   ifdef USE_GETF2_REF
    return info;
#   endif

    if (order==ColMajor) {
        for (IndexType j=0; j<min(m,n); ++j) {
//
//          Find pivot and test for singularity.
//
            IndexType jp = j + iamax(m-j, A+j*(ldA+1), IndexType(1));
            iPiv[j] = jp;
            if (A[jp+j*ldA]!=MA(0)) {
//
//              Apply the interchange to columns 1:N.
//
                if (j!=jp) {
                    swap(n, A+j, ldA, A+jp, ldA);
                }
//
//              Compute elements J+1:M of J-th column.
//
                if (j<m) {
                    scal(m-j-1, MA(1)/A[j*(ldA+1)],
                         A+j*(ldA+1)+1, IndexType(1));
                }
            } else {
                if (info==0) {
                    info = j+1;
                }
            }
//
//          Update trailing submatrix.
//
            ger(order, m-j-1, n-j-1, MA(-1),
                A+j*(ldA+1)+1, 1,
                A+j*(ldA+1)+ldA, ldA,
                A+(j+1)*(ldA+1), ldA);
        }
    }
    if (order==RowMajor) {
        assert(0);
    }

    return info;
}

} // namespace cxxlapack

#endif // CXXLAPACK_GESV_GETF2_TCC
