/*
 *   Copyright (c) 2009 Michael Lehn
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

#ifndef FLENS_BLAS_LEVEL2EXTENSIONS_MV_TCC
#define FLENS_BLAS_LEVEL2EXTENSIONS_MV_TCC 1

#include <flens/blas/closures/closures.h>
#include <flens/blas/level1/level1.h>
#include <flens/blas/level2/level2.h>
#include <flens/blas/level2extensions/level2extensions.h>
#include <flens/typedefs.h>
#include <ulmblas/cxxblas.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//-- gecrsmv
template <typename Alpha, typename MA, typename VX, typename Beta, typename VY>
typename RestrictTo<IsGeCRSMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(Transpose trans, const Alpha &alpha, const MA &A, const VX &x,
   const Beta &beta, VY &&y)
{
    const bool noTrans = (trans==NoTrans || trans==Conj);

#   ifndef NDEBUG
    if (noTrans) {
        ASSERT(x.length()==A.numCols());
    } else {
        ASSERT(x.length()==A.numRows());
    }
#   endif

    typedef typename RemoveRef<MA>::Type  MatrixA;
    typedef typename MatrixA::IndexType   IndexType;
    IndexType yLength = noTrans ? A.numRows()
                                : A.numCols();

    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(beta==Beta(0) || y.length()==yLength || y.length()==0);

    if (y.length()!=yLength) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(yLength, y.firstIndex(), Zero);
    }

//  Sparse BLAS only supports this case:
    ASSERT(x.stride()==1);
    ASSERT(y.stride()==1);

#   ifdef HAVE_CXXBLAS_GECRSMV
    cxxblas::gecrsmv(trans,
                     A.numRows(), A.numCols(),
                     alpha,
                     A.engine().values().data(),
                     A.engine().rows().data(),
                     A.engine().cols().data(),
                     x.data(),
                     beta,
                     y.data());
#   else
    ASSERT(0);
#   endif
}

//-- hecrsmv
template <typename Alpha, typename MA, typename VX, typename Beta, typename VY>
typename RestrictTo<IsHeCRSMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const Alpha &alpha, const MA &A, const VX &x, const Beta &beta, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==Beta(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

    cxxblas::hecrsmv(A.upLo(),
                     A.dim(),
                     alpha,
                     A.engine().values().data(),
                     A.engine().rows().data(),
                     A.engine().cols().data(),
                     x.data(),
                     beta,
                     y.data());
}

//-- sycrsmv
template <typename Alpha, typename MA, typename VX, typename Beta, typename VY>
typename RestrictTo<IsSyCRSMatrix<MA>::value
                 && IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
mv(const Alpha &alpha, const MA &A, const VX &x, const Beta &beta, VY &&y)
{
    ASSERT(!DEBUGCLOSURE::identical(x, y));
    ASSERT(x.length()==A.dim());
    ASSERT((beta==Beta(0)) || (y.length()==A.dim()));

    if (y.length()!=A.dim()) {
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;

        const T  Zero(0);
        y.resize(A.dim(), y.firstIndex(), Zero);
    }

    cxxblas::sycrsmv(A.upLo(),
                     A.dim(),
                     alpha,
                     A.engine().values().data(),
                     A.engine().rows().data(),
                     A.engine().cols().data(),
                     x.data(),
                     beta,
                     y.data());
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL2EXTENSIONS_MV_TCC
