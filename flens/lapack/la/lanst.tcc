/*
 *   Copyright (c) 2014, Michael Lehn
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
       DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LANST_TCC
#define FLENS_LAPACK_LA_LANST_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename VD, typename VE>
typename VD::ElementType
lanst_impl(Norm                   norm,
           const DenseVector<VD>  &d,
           const DenseVector<VE>  &e)
{
    using std::abs;
    using std::max;
    using std::sqrt;

    typedef typename VD::ElementType  T;
    typedef typename VD::IndexType    IndexType;

    const T  Zero(0), One(1);

    const Underscore<IndexType>  _;

    const IndexType n = d.length();

    T normA = 0;

    if (n==0) {
        normA = Zero;
    } else if (norm==MaximumNorm) {
//
//      Find max(abs(A(i,j))).
//
        normA = abs(d(n));
        for (IndexType i=1; i<=n-1; ++i) {
            normA = max(normA, abs(d(i)));
            normA = max(normA, abs(e(i)));
        }
    } else if (norm==OneNorm || norm==InfinityNorm) {
//
//      Find norm1(A).
//
        if (n==1) {
            normA = abs(d(1));
        } else {
            normA = max(abs(d(1))  +abs(e(1)),
                        abs(e(n-1))+abs(d(n)));
            for (IndexType i=2; i<=n-1; ++i) {
                normA = max(normA, abs(d(i))+abs(e(i))+abs(e(i-1)));
            }
        }
    } else if (norm==FrobeniusNorm) {
//
//      Find normF(A).
//
        T scale = Zero;
        T sum   = One;
        if (n>1) {
            lassq(e(_(1,n-1)), scale, sum);
            sum *= 2;
        }
        lassq(d, scale, sum);
        normA = scale*sqrt(sum);
    }
    return normA;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename VD, typename VE>
typename VD::ElementType
lanst_impl(Norm                   norm,
           const DenseVector<VD>  &d,
           const DenseVector<VE>  &e)
{
    typedef typename VD::IndexType  IndexType;

    return cxxlapack::lanst<IndexType>(getF77Char(norm),
                                       d.length(),
                                       d.data(),
                                       e.data());
}

} // namespace external

#endif

//== public interface ==========================================================

template <typename VD, typename VE>
typename RestrictTo<IsRealDenseVector<VD>::value
                 && IsRealDenseVector<VE>::value,
         typename VD::ElementType>::Type
lanst(Norm      norm,
      const VD  &d,
      const VE  &e)
{
    typedef typename VD::ElementType   T;

//
//  Test the input parameters
//
    ASSERT(d.firstIndex()==1);
    ASSERT(e.firstIndex()==1);
    ASSERT(d.length()==e.length()+1 || d.length()==0);

//
//  Call implementation
//
    T result = LAPACK_SELECT::lanst_impl(norm, d, e);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    T result_ = external::lanst_impl(norm, d, e);

    if (! isIdentical(result, result_, " result", "result_")) {
        ASSERT(0);
    }
#   endif

    return result;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LANST_TCC
