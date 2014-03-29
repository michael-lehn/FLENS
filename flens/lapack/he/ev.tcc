/*
 *   Copyright (c) 2012, Klaus Pototzky
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

/* Based on
 *
      SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
     $                  INFO )
 *
 *  -- LAPACK driver routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_HE_EV_TCC
#define FLENS_LAPACK_HE_EV_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (he)ev [real variant] -----------------------------------------------------

template <typename MA>
typename RestrictTo<IsComplex<typename MA::ElementType>::value,
         Pair<typename MA::IndexType> >::Type
ev_wsq_impl(bool                  computeV,
            HeMatrix<MA>          &A)
{
    using std::max;

    typedef typename HeMatrix<MA>::ElementType          T;
    typedef typename ComplexTrait<T>::PrimitiveType     PT;
    typedef typename HeMatrix<MA>::IndexType            IndexType;

//
//  Compute minimal workspace
//
    IndexType  n = A.dim();
    IndexType  minWork;

    if (n==0) {
        minWork = 1;
    } else {
        minWork = max(1,2*n-1);
    }

//
//  Get optimal workspace from external LAPACK
//
    T           DUMMY, WORK;
    PT          RDUMMY, RWORK;
    IndexType   LWORK = -1;
    cxxlapack::heev(computeV ? 'V' : 'N',
                    getF77Char(A.upLo()),
                    A.dim(),
                    &DUMMY,
                    A.leadingDimension(),
                    &RDUMMY,
                    &WORK,
                    LWORK,
                    &RWORK);
    return Pair<IndexType>(minWork,WORK.real());
}

//-- (he)ev [real variant] -----------------------------------------------------

template <typename MA, typename VW, typename VWORK, typename VRWORK>
typename HeMatrix<MA>::IndexType
ev_impl(bool                  computeV,
        HeMatrix<MA>          &A,
        DenseVector<VW>       &w,
        DenseVector<VWORK>    &work,
        DenseVector<VRWORK>   &rWork)
{
    using std::max;

    typedef typename HeMatrix<MA>::IndexType  IndexType;

    if (work.length()==0) {
        const auto ws = ev_wsq_impl(computeV, A);
        work.resize(ws.second, 1);
    }
    if (rWork.length()==0) {
        rWork.resize(max(1,3*A.dim()-2));
    }
    IndexType  info;
    info = cxxlapack::heev(computeV ? 'V' : 'N',
               getF77Char(A.upLo()),
                           A.dim(),
                           A.data(),
                           A.leadingDimension(),
                           w.data(),
                           work.data(),
                           work.length(),
               rWork.data());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

#ifdef USE_CXXLAPACK

//-- (he)ev [complex variant] --------------------------------------------------
template <typename MA, typename VW, typename VWORK, typename VRWORK>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsRealDenseVector<VW>::value
                 && IsComplexDenseVector<VWORK>::value
                 && IsRealDenseVector<VRWORK>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
ev(bool     computeV,
   MA       &&A,
   VW       &&w,
   VWORK    &&work,
   VRWORK   &&rWork)
{
    LAPACK_DEBUG_OUT("(he)ev [complex]");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type      MatrixA;
    typedef typename MatrixA::IndexType       IndexType;

    const IndexType n = A.dim();

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(work.firstIndex()==1);
    ASSERT(rWork.firstIndex()==1);

    ASSERT(w.firstIndex()==1);
    ASSERT(w.length()==0 || w.length()==n);

#   endif

//
//  Resize output arguments if they are empty and needed
//
    if (w.length()==0) {
        w.resize(n, 1);
    }

//
//  Call external implementation
//
    IndexType result = external::ev_impl(computeV, A, w, work, rWork);
    return result;
}

//-- (he)ev [real variant with temporary workspace] ----------------------------

template <typename MA, typename VW>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsRealDenseVector<VW>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
ev(bool     computeV,
   MA       &&A,
   VW       &&w)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type::Vector        WorkVector;
    typedef typename RemoveRef<MA>::Type::ElementType   T;
    typedef typename ComplexTrait<T>::PrimitiveType     PT;
    typedef DenseVector<Array<PT> >                     RealWorkVector;

    WorkVector      work;
    RealWorkVector  rWork;

    return ev(computeV, A, w, work, rWork);
}

#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_HE_EV_TCC
