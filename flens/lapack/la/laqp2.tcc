/*
 *   Copyright (c) 2012, Michael Lehn
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
       SUBROUTINE DLAQP2( M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2,
      $                   WORK )
 *
 *  -- LAPACK auxiliary routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_LA_LAQP2_TCC
#define FLENS_LAPACK_LA_LAQP2_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//
//  Real variant
//
template <typename MA, typename JPIV, typename VTAU,
          typename VN1, typename VN2, typename VWORK>
typename RestrictTo<IsReal<typename MA::ElementType>::value,
         void>::Type
laqp2_impl(typename GeMatrix<MA>::IndexType  offset,
           GeMatrix<MA>                      &A,
           DenseVector<JPIV>                 &jPiv,
           DenseVector<VTAU>                 &tau,
           DenseVector<VN1>                  &vn1,
           DenseVector<VN2>                  &vn2,
           DenseVector<VWORK>                &work)
{
    using std::abs;
    using std::max;
    using std::min;
    using flens::pow;
    using std::sqrt;
    using std::swap;

    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const Underscore<IndexType> _;

    const IndexType m  = A.numRows();
    const IndexType n  = A.numCols();

    const IndexType mn = min(m-offset, n);

    const ElementType Zero(0), One(1);
    const ElementType tol3z = sqrt(lamch<ElementType>(Eps));
//
//  Compute factorization.
//
    for (IndexType i=1; i<=mn; ++i) {

        IndexType offpi = offset + i;
//
//      Determine ith pivot column and swap if necessary.
//
        IndexType pvt = (i-1) + blas::iamax(vn1(_(i,n)));

        if (pvt!=i) {
            blas::swap(A(_,pvt), A(_,i));
            swap(jPiv(pvt), jPiv(i));
            vn1(pvt) = vn1(i);
            vn2(pvt) = vn2(i);
        }
//
//      Generate elementary reflector H(i).
//
        if (offpi<m) {
            larfg(m-offpi+1, A(offpi,i), A(_(offpi+1,m),i), tau(i));
        } else {
            larfg(1, A(m,i), A(_(m+1,m),i), tau(i));
        }

        if (i<=n) {
//
//          Apply H(i)**T to A(offset+i:m,i+1:n) from the left.
//
            const ElementType Aii = A(offpi,i);
            A(offpi,i) = One;
            larf(Left, A(_(offpi,m),i), cxxblas::conjugate(tau(i)), A(_(offpi,m),_(i+1,n)),
                 work(_(1,n)));
            A(offpi,i) = Aii;
        }
//
//      Update partial column norms.
//
        for (IndexType j=i+1; j<=n; ++j) {
            if (vn1(j)!=Zero) {
//
//              NOTE: The following 4 lines follow from the analysis in
//              Lapack Working Note 176.
//
                ElementType tmp = One - pow(abs(A(offpi,j))/vn1(j), 2);
                tmp = max(tmp, Zero);
                ElementType tmp2 = tmp*pow(vn1(j)/vn2(j), 2);
                if (tmp2<=tol3z) {
                    if (offpi<m) {
                        vn1(j) = blas::nrm2(A(_(offpi+1,m),j));
                        vn2(j) = vn1(j);
                    } else {
                        vn1(j) = Zero;
                        vn2(j) = Zero;
                    }
                } else {
                    vn1(j) *= sqrt(tmp);
                }
            }
        }

    }
}

//
//  Complex variant
//
template <typename MA, typename JPIV, typename VTAU,
          typename VN1, typename VN2, typename VWORK>
typename RestrictTo<IsComplex<typename MA::ElementType>::value,
         void>::Type
laqp2_impl(typename GeMatrix<MA>::IndexType  offset,
           GeMatrix<MA>                      &A,
           DenseVector<JPIV>                 &jPiv,
           DenseVector<VTAU>                 &tau,
           DenseVector<VN1>                  &vn1,
           DenseVector<VN2>                  &vn2,
           DenseVector<VWORK>                &work)
{
    using std::abs;
    using std::conj;
    using std::max;
    using std::min;
    using flens::pow;
    using std::sqrt;
    using std::swap;

    typedef typename GeMatrix<MA>::ElementType                 ElementType;
    typedef typename ComplexTrait<ElementType>::PrimitiveType  PrimitiveType;
    typedef typename GeMatrix<MA>::IndexType                   IndexType;

    const Underscore<IndexType> _;

    const IndexType m  = A.numRows();
    const IndexType n  = A.numCols();

    const IndexType mn = min(m-offset, n);

    const PrimitiveType  Zero(0), One(1);
    const ElementType    CZero(0), COne(1);
    const PrimitiveType  tol3z = sqrt(lamch<PrimitiveType>(Eps));
//
//  Compute factorization.
//
    for (IndexType i=1; i<=mn; ++i) {

        IndexType offpi = offset + i;
//
//      Determine ith pivot column and swap if necessary.
//
        IndexType pvt = (i-1) + blas::iamax(vn1(_(i,n)));

        if (pvt!=i) {
            blas::swap(A(_,pvt), A(_,i));
            swap(jPiv(pvt), jPiv(i));
            vn1(pvt) = vn1(i);
            vn2(pvt) = vn2(i);
        }
//
//      Generate elementary reflector H(i).
//
        if (offpi<m) {
            larfg(m-offpi+1, A(offpi,i), A(_(offpi+1,m),i), tau(i));
        } else {
            larfg(1, A(m,i), A(_(m+1,m),i), tau(i));
        }

        if (i<=n) {
//
//          Apply H(i)**T to A(offset+i:m,i+1:n) from the left.
//
            const ElementType Aii = A(offpi,i);
            A(offpi,i) = One;
            larf(Left, A(_(offpi,m),i), conj(tau(i)), A(_(offpi,m),_(i+1,n)),
                 work(_(1,n)));
            A(offpi,i) = Aii;
        }
//
//      Update partial column norms.
//
        for (IndexType j=i+1; j<=n; ++j) {
            if (vn1(j)!=Zero) {
//
//              NOTE: The following 4 lines follow from the analysis in
//              Lapack Working Note 176.
//
                PrimitiveType tmp = One - pow(abs(A(offpi,j))/vn1(j), 2);
                tmp = max(tmp, Zero);
                PrimitiveType tmp2 = tmp*pow(vn1(j)/vn2(j), 2);
                if (tmp2<=tol3z) {
                    if (offpi<m) {
                        vn1(j) = blas::nrm2(A(_(offpi+1,m),j));
                        vn2(j) = vn1(j);
                    } else {
                        vn1(j) = Zero;
                        vn2(j) = Zero;
                    }
                } else {
                    vn1(j) *= sqrt(tmp);
                }
            }
        }

    }
}


} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename JPIV, typename VTAU,
          typename VN1, typename VN2, typename VWORK>
void
laqp2_impl(typename GeMatrix<MA>::IndexType  offset,
           GeMatrix<MA>                      &A,
           DenseVector<JPIV>                 &jPiv,
           DenseVector<VTAU>                 &tau,
           DenseVector<VN1>                  &vn1,
           DenseVector<VN2>                  &vn2,
           DenseVector<VWORK>                &work)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    cxxlapack::laqp2<IndexType>(A.numRows(),
                                A.numCols(),
                                offset,
                                A.data(),
                                A.leadingDimension(),
                                jPiv.data(),
                                tau.data(),
                                vn1.data(),
                                vn2.data(),
                                work.data());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
//
//  Real and complex variant
//
template <typename IndexType, typename MA, typename JPIV, typename VTAU,
          typename VN1, typename VN2, typename VWORK>
    typename RestrictTo<((IsRealGeMatrix<MA>::value &&
                          IsRealDenseVector<VTAU>::value &&
                          IsRealDenseVector<VWORK>::value)
                      || (IsComplexGeMatrix<MA>::value &&
                          IsComplexDenseVector<VTAU>::value &&
                          IsComplexDenseVector<VWORK>::value))
                      && IsIntegerDenseVector<JPIV>::value
                      && IsRealDenseVector<VN1>::value
                      && IsRealDenseVector<VN2>::value,
             void>::Type
laqp2(IndexType  offset,
      MA         &&A,
      JPIV       &&jPiv,
      VTAU       &&tau,
      VN1        &&vn1,
      VN2        &&vn2,
      VWORK      &&work)
{
    using std::min;

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MA>::Type        MatrixA;
    typedef typename RemoveRef<JPIV>::Type      VectorJPiv;
    typedef typename RemoveRef<VTAU>::Type      VectorTau;
    typedef typename RemoveRef<VN1>::Type       VectorVN1;
    typedef typename RemoveRef<VN2>::Type       VectorVN2;
    typedef typename RemoveRef<VWORK>::Type     VectorWork;
#   endif


#   ifndef NDEBUG
//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(jPiv.firstIndex()==1);
    ASSERT(tau.firstIndex()==1);
    ASSERT(vn1.firstIndex()==1);
    ASSERT(vn2.firstIndex()==1);
    ASSERT(work.firstIndex()==1);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    // Lehn: I think there is a bug in the DLAQP2 documentation.  This should
    //       be the length of tau.
    const IndexType k = min(m-offset, n);

    ASSERT(jPiv.length()==n);
    ASSERT(tau.length()==k);
    ASSERT(vn1.length()==n);
    ASSERT(vn1.length()==n);
    ASSERT(work.length()==n);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename MatrixA::NoView     A_org      = A;
    typename VectorJPiv::NoView  jPiv_org   = jPiv;
    typename VectorTau::NoView   tau_org    = tau;
    typename VectorVN1::NoView   vn1_org    = vn1;
    typename VectorVN2::NoView   vn2_org    = vn2;
    typename VectorWork::NoView  work_org   = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::laqp2_impl(offset, A, jPiv, tau, vn1, vn2, work);

#   ifdef CHECK_CXXLAPACK
//
//  Restore output arguments
//
    typename MatrixA::NoView     A_generic    = A;
    typename VectorJPiv::NoView  jPiv_generic = jPiv;
    typename VectorTau::NoView   tau_generic  = tau;
    typename VectorVN1::NoView   vn1_generic  = vn1;
    typename VectorVN2::NoView   vn2_generic  = vn2;
    typename VectorWork::NoView  work_generic = work;

    A    = A_org;
    jPiv = jPiv_org;
    tau  = tau_org;
    vn1  = vn1_org;
    vn2  = vn2_org;
    work = work_org;

//
//  Compare results
//
    external::laqp2_impl(offset, A, jPiv, tau, vn1, vn2, work);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(jPiv_generic, jPiv, "jPiv_generic", "jPiv")) {
        std::cerr << "CXXLAPACK: jPiv_generic = " << jPiv_generic << std::endl;
        std::cerr << "F77LAPACK: jPiv = " << jPiv << std::endl;
        failed = true;
    }

    if (! isIdentical(tau_generic, tau, "tau_generic", "tau")) {
        std::cerr << "CXXLAPACK: tau_generic = " << tau_generic << std::endl;
        std::cerr << "F77LAPACK: tau = " << tau << std::endl;
        failed = true;
    }

    if (! isIdentical(vn1_generic, vn1, "vn1_generic", "vn1")) {
        std::cerr << "CXXLAPACK: vn1_generic = " << vn1_generic << std::endl;
        std::cerr << "F77LAPACK: vn1 = " << vn1 << std::endl;
        failed = true;
    }

    if (! isIdentical(vn2_generic, vn2, "vn2_generic", "vn2")) {
        std::cerr << "CXXLAPACK: vn2_generic = " << vn2_generic << std::endl;
        std::cerr << "F77LAPACK: vn2 = " << vn2 << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "offset = " << offset << std::endl;
        std::cerr << "A = " << A << std::endl;
        std::cerr << "jPiv = " << jPiv << std::endl;
        std::cerr << "tau = " << tau << std::endl;
        std::cerr << "vn1 = " << vn1 << std::endl;
        std::cerr << "vn2 = " << vn2 << std::endl;
        std::cerr << "work =  " << work << std::endl;
        ASSERT(0);
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAQP2_TCC
