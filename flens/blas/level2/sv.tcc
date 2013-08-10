/*
 *   Copyright (c) 2010, Michael Lehn
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

#ifndef FLENS_BLAS_LEVEL2_SV_TCC
#define FLENS_BLAS_LEVEL2_SV_TCC

#include <flens/blas/closures/closures.h>
#include <flens/blas/level2/level2.h>
#include <flens/typedefs.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif


namespace flens { namespace blas {

//-- tbsv
template <typename MA, typename VX>
typename RestrictTo<IsTbMatrix<MA>::value
                 && IsDenseVector<VX>::value,
         void>::Type
sv(Transpose trans, const MA &A, VX &&x)
{
    ASSERT(x.length()==A.dim());
#   ifdef HAVE_CXXBLAS_TBSV
    cxxblas::tbsv(A.order(), A.upLo(),
                  trans, A.diag(),
                  A.dim(), A.numOffDiags(),
                  A.data(), A.leadingDimension(),
                  x.data(), x.stride());
#   else
    ASSERT(0);
#   endif
}

//-- trsv
template <typename MA, typename VX>
typename RestrictTo<IsTrMatrix<MA>::value
                 && IsDenseVector<VX>::value,
         void>::Type
sv(Transpose trans, const MA &A, VX &&x)
{

    ASSERT(x.length()==A.dim());
#   ifdef HAVE_CXXBLAS_TRSV
    cxxblas::trsv(A.order(), A.upLo(),
                  trans, A.diag(),
                  A.dim(),
                  A.data(), A.leadingDimension(),
                  x.data(), x.stride());
#   else
    ASSERT(0);
#   endif
}

//-- tpsv
template <typename MA, typename VX>
typename RestrictTo<IsTpMatrix<MA>::value
                 && IsDenseVector<VX>::value,
         void>::Type
sv(Transpose trans, const MA &A, VX &&x)
{
    ASSERT(x.length()==A.dim());
#   ifdef HAVE_CXXBLAS_TPSV
    cxxblas::tpsv(A.order(), A.upLo(),
                  trans, A.diag(),
                  A.dim(),
                  A.data(),
                  x.data(), x.stride());
#   else
    ASSERT(0);
#   endif
}

//-- trccssv
template <typename ALPHA, typename MA, typename VX, typename VY>
typename RestrictTo<IsTrCCSMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
sv(Transpose trans, const ALPHA &alpha, const MA &A, const VX &x, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));

    if (y.length()==0) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

//  Sparse BLAS only supports this case:
    ASSERT(x.stride()==1);
    ASSERT(y.stride()==1);
    
    ASSERT(x.length()==A.dim());
    ASSERT(x.length()==y.length());
    
    
    
#   ifdef HAVE_CXXBLAS_TRCCSSV
    cxxblas::trccssv(A.upLo(), trans, 
                     A.dim(),
                     alpha,
                     A.engine().values().data(),
                     A.engine().rows().data(),
                     A.engine().cols().data(),
                     x.data(),
                     y.data());
#   else
    ASSERT(0);
#   endif
    
}

//-- trcrssv
template <typename ALPHA, typename MA, typename VX, typename VY>
typename RestrictTo<IsTrCRSMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
sv(Transpose trans, const ALPHA &alpha, const MA &A, const VX &x, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));

    if (y.length()==0) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

//  Sparse BLAS only supports this case:
    ASSERT(x.stride()==1);
    ASSERT(y.stride()==1);
    
    ASSERT(x.length()==A.dim());
    ASSERT(x.length()==y.length());
    
    
    
#   ifdef HAVE_CXXBLAS_TRCRSSV
    cxxblas::trcrssv(A.upLo(), trans, 
                     A.dim(),
                     alpha,
                     A.engine().values().data(),
                     A.engine().rows().data(),
                     A.engine().cols().data(),
                     x.data(),
                     y.data());
#   else
    ASSERT(0);
#   endif
    
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL3_SV_TCC
