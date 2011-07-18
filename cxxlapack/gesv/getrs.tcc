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

#ifndef CXXLAPACK_GESV_GETRS_TCC
#define CXXLAPACK_GESV_GETRS_TCC 1

#include <cxxblas/cxxblas.h>
#include <cxxlapack/cxxlapack.h>

namespace cxxlapack {

using cxxblas::ger;
using cxxblas::iamax;
using cxxblas::scal;
using cxxblas::swap;

template <typename IndexType, typename MA, typename MB>
IndexType
getrs(StorageOrder order, Transpose transA,
      IndexType n, IndexType nRhs,
      const MA *A, IndexType ldA,
      const IndexType *iPiv,
      MB *B, IndexType ldB)
{
    IndexType info = 0;

    if (n<IndexType(0)) {
        info = -2;
        // xerbla("getrs", info);
        return info;
    }
    if (nRhs<IndexType(0)) {
        info = -3;
        // xerbla("getrs", info);
        return info;
    }
    if (ldA<max(IndexType(1), n)) {
        info = -5;
        // xerbla("getrs", info);
        return info;
    }
    if (ldB<max(IndexType(1), n)) {
        info = -8;
        // xerbla("getrs", info);
        return info;
    }

//
//     Quick return if possible
//
    if ((n==0) || (nRhs==0)) {
        return info;
    }

#   ifdef USE_GETRS_REF
    return info;
#   endif

    if ((transA==NoTrans) || (transA==Conj)) {
//
//      Solve A * X = B.
//
//      Apply row interchanges to the right hand sides.
//
        laswp(order, nRhs, B, ldB, IndexType(0), n-1, iPiv, IndexType(1));
//
//      Solve L*X = B, overwriting B with X.
//
        trsm(order, Left, Lower, NoTrans, Unit,
             n, nRhs, MA(1), A, ldA, B, ldB);
//
//      Solve U*X = B, overwriting B with X.
//
        trsm(order, Left, Upper, NoTrans, NonUnit,
             n, nRhs, MA(1), A, ldA, B, ldB);
    } else {
//
//      Solve A' * X = B.
//
//      Solve U'*X = B, overwriting B with X.
//
        trsm(order, Left, Upper, transA, NonUnit,
             n, nRhs, MA(1), A, ldA, B, ldB);
//
//      Solve L'*X = B, overwriting B with X.
//
        trsm(order, Left, Lower, transA, Unit,
             n, nRhs, MA(1), A, ldA, B, ldB);
//
//      Apply row interchanges to the solution vectors.
//
        laswp(order, nRhs, B, ldB, IndexType(0), n-1, iPiv, IndexType(-1));
    }
    return info;
}

} // namespace cxxlapack

#endif // CXXLAPACK_GESV_GETRS_TCC
