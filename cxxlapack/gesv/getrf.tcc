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

#ifndef CXXLAPACK_GESV_GETRF_TCC
#define CXXLAPACK_GESV_GETRF_TCC 1

#include <algorithm>
#include <cxxblas/cxxblas.h>
#include <cxxlapack/cxxlapack.h>

namespace cxxlapack {

using cxxblas::gemm;
using cxxblas::iamax;
using cxxblas::scal;
using cxxblas::swap;

template <typename IndexType, typename MA>
IndexType
getrf(StorageOrder order,
      IndexType m, IndexType n, MA *A, IndexType ldA,
      IndexType *iPiv)
{
//
//  Test the input parameters.
//
    IndexType info = 0;
    if (m<0) {
        info = -1;
        // xerbla("getrf", info);
        return info;
    }
    if (n<0) {
        info = -2;
        // xerbla("getrf", info);
        return info;
    }
    if (ldA<max(1,m)) {
        info = -4;
        // xerbla("getrf", info);
        return info;
    }
//
//  Quick return if possible
//
    if ((m==0) || (n==0)) {
        return info;
    }
    
#   ifdef USE_GETRF_REF
    return info;
#   endif

//
//  Determine the block size for this environment.
//
    IndexType nb = 2; // laenv( 1, 'DGETRF', ' ', M, N, -1, -1 );
    if ((nb<=1) || (nb>=min(m,n))) {
//
//      Use unblocked code.
//
        info = getf2(order, m, n, A, ldA, iPiv);
    } else {
//
//      Use blocked code.
//
        for (IndexType j=1; j<=min(m,n); j+=bs) {
            IndexType jb = min(min(m,n)-j, bs);
//
//          Factor diagonal and subdiagonal blocks and test for exact
//          singularity.
//
            Range rows(j, m);
            Range cols(j, j+jb);
            IndexType _info = getf2(A(rows,cols), iPiv(rows));
//
//          Adjust INFO and the pivot indices.
//
            if ((info==0) && (_info>0)) {
                info = _info + j;
            }
            for (IndexType i=j; i<min(m,j+jb); ++i) {
                iPiv[i] += j;
            }
//
//          Apply interchanges to columns 0:J.
//
            laswp(order, j, A, ldA, j, j+jb-1, iPiv, IndexType(1));

            if (j+jb<n) {
//
//              Apply interchanges to columns J+JB:N.
//
                laswp(order, n-j-jb+1, A+(j+jb)*ldA, ldA,
                             j, j+jb-1, iPiv, IndexType(1));
//
//              Compute block row of U.
//
                trsm(order, Left, Lower, NoTrans, Unit,
                     jb, n-j-jb,
                     MA(1),
                     A + j*(ldA+1), ldA,
                     A + j+(j+jb)*ldA, ldA);

                if (j+jb<m) {
//
//                  Update trailing submatrix.
//
                    gemm(order, NoTrans, NoTrans,
                         m-j-jb, n-j-jb, jb,
                         MA(-1),
                         A+j*(ldA+1)+jb, ldA,
                         A+j+(j+jb)*ldA, ldA,
                         MA(1),
                         A+(j+jb)*(ldA+1), ldA);
                }
            }
        }
    }
    return info;
}

} // namespace cxxlapack

#endif // CXXLAPACK_GESV_GETRF_TCC
