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

#ifndef CXXLAPACK_GESV_GEEQU_TCC
#define CXXLAPACK_GESV_GEEQU_TCC 1

#include <cxxblas/cxxblas.h>
#include <cxxlapack/cxxlapack.h>

namespace cxxlapack {

template <typename IndexType, typename T>
IndexType
geequ(StorageOrder order,
      IndexType m, IndexType n,
      const T *A, IndexType ldA,
      T *r, T *c,
      T &rowCond, T &colCond,
      T &AMax)
{
    IndexType info = 0;

    if (m<0) {
        info = -1;
    } else if (n<0) {
        info = -2;
    } else if (ldA<max(1,m)) {
        info = -4;
    }
    if (info!=0) {
        // xerbla
        return info;
    }
//
//  Quick return if possible
//
    if ((m==0) || (n==0)) {
        rowCond = 1;
        colCond = 1;
        AMax = 0;
        return info;
    }
#   ifdef SHOW_TODO
    fprintf(stderr, "TODO: GEEQU\n");
#   endif
    return info;
//
//  Get machine constants.
//
    const T smallNum = lamch<T>(SafeMin);
    const T bigNum = T(1) / smallNum;
//
//  Compute row scale factors.
//
    for (IndexType i=0; i<m; ++i) {
        r[i] = T(0);
    }
//
//  Find the maximum element in each row.
//
    for (IndexType j=0; j<n; ++j) {
        for (IndexType i=0; i<m; ++i) {
            r[i] = max(r[i], abs(A[i+ldA*j]));
        }
    }
//
//  Find the maximum and minimum scale factors.
//
    T rcMin = bigNum;
    T rcMax = 0;

    for (IndexType i=0; i<m; ++i) {
        rcMax = max(rcMax, r[i]);
        rcMin = min(rcMin, r[i]);
    }
    AMax = rcMax;

    if (rcMin==T(0)) {
//
//      Find the first zero scale factor and return an error code.
//
        for (IndexType i=0; i<m; ++i) {
            if (r[i]==T(0)) {
                return i+1;
            }
        }
    } else {
//
//      Invert the scale factors.
//
        for (IndexType i=0; i<m; ++i) {
            r[i] = T(1) / min(max(r[i],smallNum), bigNum);
        }
//
//      Compute ROWCND = min(R(I)) / max(R(I))
//
        rowCond = max(rcMin, smallNum) / min(rcMax, bigNum);
    }
//
//    Compute column scale factors
//
    for (IndexType j=0; j<n; ++j) {
        c[j] = 0;
    }
//
//  Find the maximum element in each column,
//  assuming the row scaling computed above.
//
    for (IndexType j=0; j<n; ++j) {
        for (IndexType i=0; i<m; ++i) {
            c[j] = max(c[j],abs(A[i+j*ldA]))*r[i];
        }
    }
//
//  Find the maximum and minimum scale factors.
//
    rcMin = bigNum;
    rcMax = 0;
    for (IndexType j=0; j<n; ++j) {
        rcMin = min(rcMin, c[j]);
        rcMax = max(rcMax, c[j]);
    }

    if (rcMin==T(0)) {
//
//      Find the first zero scale factor and return an error code.
//
        for (IndexType j=0; j<n; ++j) {
            if (c[j]==T(0)) {
                return m+j+1;
            }
        }
    } else {
//
//      Invert the scale factors.
//
        for (IndexType j=0; j<n; ++j) {
            c[j] = T(1) / min(max(c[j],smallNum), bigNum);
        }
//
//      Compute COLCND = min(C(J)) / max(C(J))
//
        colCond = max(rcMin, smallNum) / min(rcMax, bigNum);
    }
    return info;
}

} // namespace cxxlapack

#endif // CXXLAPACK_GESV_GEEQU_TCC
