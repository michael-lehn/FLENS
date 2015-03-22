/*
 *   Copyright (c) 2011, Michael Lehn
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
      SUBROUTINE DLAQR1( N, H, LDH, SR1, SI1, SR2, SI2, V )
      SUBROUTINE ZLAQR1( N, H, LDH, S1, S2, V )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *     Univ. of Tennessee, Univ. of California Berkeley,
 *     Univ. of Colorado Denver and NAG Ltd..
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LAQR1_TCC
#define FLENS_LAPACK_LA_LAQR1_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//
//  Real variant
//
template <typename MH, typename T, typename VV>
void
laqr1_impl(const GeMatrix<MH>        &H,
           const T                   &sr1,
           const T                   &si1,
           const T                   &sr2,
           const T                   &si2,
           DenseVector<VV>           &v)
{
    using std::abs;

    typedef typename GeMatrix<MH>::IndexType    IndexType;

    const IndexType n   = H.numRows();
    const T         Zero(0);

    if (n==2) {
        const T s = abs(H(1,1)-sr2) + abs(si2) + abs(H(2,1));
        if (s==Zero) {
            v(1) = Zero;
            v(2) = Zero;
        } else {
            const T H21s = H(2,1)/s;
            v(1) = H21s*H(1,2) + (H(1,1)-sr1)*((H(1,1)-sr2)/s) - si1*(si2/s);
            v(2) = H21s*(H(1,1) + H(2,2)-sr1-sr2);
        }
    } else {
        const T s = abs(H(1,1)-sr2) + abs(si2) + abs(H(2,1)) + abs(H(3,1));
        if (s==Zero) {
            v(1) = Zero;
            v(2) = Zero;
            v(3) = Zero;
        } else {
            const T H21s = H(2,1) / s;
            const T H31s = H(3,1) / s;
            v(1) = (H(1,1)-sr1)*((H(1,1)-sr2)/s) - si1*(si2/s)
                    + H(1,2)*H21s + H(1,3)*H31s;
            v(2) = H21s*(H(1,1) + H(2,2)-sr1-sr2) + H(2,3 )*H31s;
            v(3) = H31s*(H(1,1) + H(3,3)-sr1-sr2) + H21s*H(3,2);
        }
    }
}

//
//  Complex variant
//
template <typename MH, typename T, typename VV>
void
laqr1_impl(const GeMatrix<MH>        &H,
           const T                   &s1,
           const T                   &s2,
           DenseVector<VV>           &v)
{
    typedef typename GeMatrix<MH>::IndexType    IndexType;

    const IndexType n   = H.numRows();
    const T         Zero(0);

    if (n==2) {
        const T s = abs1(H(1,1)-s2) + abs1(H(2,1));
        if (s==Zero) {
            v(1) = Zero;
            v(2) = Zero;
        } else {
            const T H21s = H(2,1)/s;
            v(1) = H21s*H(1,2) + (H(1,1)-s1)*((H(1,1)-s2)/s);
            v(2) = H21s*(H(1,1) + H(2,2) - s1 - s2);
        }
    } else {
        const T s = abs1(H(1,1)-s2) + abs1(H(2,1)) + abs1(H(3,1));
        if (s==Zero) {
            v(1) = Zero;
            v(2) = Zero;
            v(3) = Zero;
        } else {
            const T H21s = H(2,1) / s;
            const T H31s = H(3,1) / s;
            v(1) = (H(1,1)-s1)*((H(1,1)-s2)/s) + H(1,2)*H21s + H(1,3)*H31s;
            v(2) = H21s*(H(1,1) + H(2,2)-s1-s2) + H(2,3)*H31s;
            v(3) = H31s*(H(1,1) + H(3,3)-s1-s2) + H21s*H(3,2);
        }
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//
//  Real variant
//
template <typename MH, typename T, typename VV>
void
laqr1_impl(const GeMatrix<MH>        &H,
           const T                   &sr1,
           const T                   &si1,
           const T                   &sr2,
           const T                   &si2,
           DenseVector<VV>           &v)
{
    typedef typename GeMatrix<MH>::IndexType  IndexType;

    cxxlapack::laqr1<IndexType>(H.numRows(),
                                H.data(),
                                H.leadingDimension(),
                                sr1,
                                si1,
                                sr2,
                                si2,
                                v.data());
}

//
//  Complex variant
//
template <typename MH, typename T, typename VV>
void
laqr1_impl(const GeMatrix<MH>        &H,
           const T                   &s1,
           const T                   &s2,
           DenseVector<VV>           &v)
{
    typedef typename GeMatrix<MH>::IndexType  IndexType;

    cxxlapack::laqr1<IndexType>(H.numRows(),
                                H.data(),
                                H.leadingDimension(),
                                s1,
                                s2,
                                v.data());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
//
//  Real variant
//
template <typename MH, typename T, typename VV>
typename RestrictTo<IsRealGeMatrix<MH>::value
                 && IsRealDenseVector<VV>::value,
         void>::Type
laqr1(const MH  &H,
      const T   &sr1,
      const T   &si1,
      const T   &sr2,
      const T   &si2,
      VV        &&v)
{
    LAPACK_DEBUG_OUT("laqr1");

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<VV>::Type   VectorV;
#   endif

//
//  Test the input parameters
//
#   ifndef NDEBUG
    typedef typename RemoveRef<MH>::Type   MatrixH;
    typedef typename MatrixH::IndexType    IndexType;

    ASSERT(H.firstRow()==1);
    ASSERT(H.firstCol()==1);
    ASSERT(H.numRows()==H.numCols());

    const IndexType n = H.numRows();

    ASSERT(n==2 || n==3);

    ASSERT(v.length()==n);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename VectorV::NoView    v_org  = v;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::laqr1_impl(H, sr1, si1, sr2, si2, v);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename VectorV::NoView    v_generic  = v;
//
//  restore output arguments
//
    v = v_org;
//
//  Compare results
//
    external::laqr1_impl(H, sr1, si1, sr2, si2, v);

    bool failed = false;

    if (! isIdentical(v_generic, v, "v_generic", "v")) {
        std::cerr << "CXXLAPACK: v_generic = " << v_generic << std::endl;
        std::cerr << "F77LAPACK: v = " << v << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: laqr1.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: laqr1.tcc" << std::endl;
    }
#   endif
}

//
//  Complex variant
//
template <typename MH, typename T, typename VV>
typename RestrictTo<IsComplexGeMatrix<MH>::value
                 && IsComplexDenseVector<VV>::value,
         void>::Type
laqr1(const MH  &H,
      const T   &s1,
      const T   &s2,
      VV        &&v)
{
    LAPACK_DEBUG_OUT("laqr1 (complex)");

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<VV>::Type   VectorV;
#   endif
//
//  Test the input parameters
//
#   ifndef NDEBUG
    typedef typename RemoveRef<MH>::Type   MatrixH;
    typedef typename MatrixH::IndexType    IndexType;

    ASSERT(H.firstRow()==1);
    ASSERT(H.firstCol()==1);
    ASSERT(H.numRows()==H.numCols());

    const IndexType n = H.numRows();

    ASSERT(n==2 || n==3);

    ASSERT(v.length()==n);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename VectorV::NoView    v_org  = v;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::laqr1_impl(H, s1, s2, v);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename VectorV::NoView    v_generic  = v;
//
//  restore output arguments
//
    v = v_org;
//
//  Compare results
//
    external::laqr1_impl(H, s1, s2, v);

    bool failed = false;
    if (! isIdentical(v_generic, v, "v_generic", "v")) {
        std::cerr << "CXXLAPACK: v_generic = " << v_generic << std::endl;
        std::cerr << "F77LAPACK: v = " << v << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: laqr1.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: laqr1.tcc" << std::endl;
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAQR1_TCC
