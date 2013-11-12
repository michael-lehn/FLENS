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
       SUBROUTINE DLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1,
      $                   VN2, AUXV, F, LDF )
 *
 *  -- LAPACK auxiliary routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_LA_LAQPS_TCC
#define FLENS_LAPACK_LA_LAQPS_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//
//  Real variant
//
template <typename MA, typename JPIV, typename VTAU,
          typename VN1, typename VN2, typename VAUX,
          typename MF>
typename RestrictTo<IsReal<typename MA::ElementType>::value,
         void>::Type
laqps_impl(typename GeMatrix<MA>::IndexType  offset,
           typename GeMatrix<MA>::IndexType  nb,
           typename GeMatrix<MA>::IndexType  &kb,
           GeMatrix<MA>                      &A,
           DenseVector<JPIV>                 &jPiv,
           DenseVector<VTAU>                 &tau,
           DenseVector<VN1>                  &vn1,
           DenseVector<VN2>                  &vn2,
           DenseVector<VAUX>                 &aux,
           GeMatrix<MF>                      &F)
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

    const IndexType m     = A.numRows();
    const IndexType n     = A.numCols();

    const IndexType lastRk = min(m, n+offset);

    IndexType lasticc = 0;
    IndexType k = 0;

    const ElementType Zero(0), One(1);
    const ElementType tol3z = sqrt(lamch<ElementType>(Eps));
//
//  Beginning of while loop.
//
    while (k<nb && lasticc==0) {
        ++k;

        const IndexType rk = offset + k;
//
//      Determine ith pivot column and swap if necessary
//
        IndexType pvt = (k-1) + blas::iamax(vn1(_(k,n)));
        if (pvt!=k) {
            blas::swap(A(_,pvt), A(_,k));
            blas::swap(F(pvt,_(1,k-1)), F(k,_(1,k-1)));
            swap(jPiv(pvt),jPiv(k));
            vn1(pvt) = vn1(k);
            vn2(pvt) = vn2(k);
        }
//
//      Apply previous Householder reflectors to column K:
//      A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)**T.
//
        if (k>1) {
            blas::mv(NoTrans, -One, A(_(rk,m),_(1,k-1)), F(k,_(1,k-1)),
                     One, A(_(rk,m),k));
        }
//
//      Generate elementary reflector H(k).
//
        if (rk<m) {
            larfg(m-rk+1, A(rk,k), A(_(rk+1,m),k), tau(k));
        } else {
            larfg(1, A(m,k), A(_(m+1,m),k), tau(k));
        }

        const ElementType Akk = A(rk,k);
        A(rk,k) = One;
//
//      Compute Kth column of F:
//
//      Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**T*A(RK:M,K).
//
        if (k<n) {
            blas::mv(Trans, tau(k), A(_(rk,m),_(k+1,n)), A(_(rk,m),k),
                     Zero, F(_(k+1,n),k));
        }
//
//      Padding F(1:K,K) with zeros.
//
        F(_(1,k),k) = Zero;
//
//      Incremental updating of F:
//      F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)**T
//                   *A(RK:M,K).
//
        if (k>1) {
            blas::mv(Trans, -tau(k), A(_(rk,m),_(1,k-1)), A(_(rk,m),k),
                     Zero, aux(_(1,k-1)));
            blas::mv(NoTrans, One, F(_,_(1,k-1)), aux(_(1,k-1)),
                     One, F(_,k));
        }
//
//      Update the current row of A:
//      A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)**T.
//
        if (k<n) {
            blas::mv(NoTrans, -One, F(_(k+1,n),_(1,k)), A(rk,_(1,k)),
                     One, A(rk,_(k+1,n)));
        }
//
//      Update partial column norms.
//
        if (rk<lastRk) {
            for (IndexType j=k+1; j<=n; ++j) {
                if (vn1(j)!=Zero) {
//
//                  NOTE: The following 4 lines follow from the analysis in
//                  Lapack Working Note 176.
//
                    ElementType tmp = abs(A(rk,j)) / vn1(j);
                    tmp = max(Zero, (One+tmp)*(One-tmp));
                    ElementType tmp2 = tmp*pow(vn1(j)/vn2(j),2);
                    if (tmp2<=tol3z) {
                        vn2(j) = ElementType(lasticc);
                        lasticc = j;
                    } else {
                        vn1(j) *= sqrt(tmp);
                    }
                }
            }
        }

        A(rk,k) = Akk;
//
//      End of while loop.
//
    }
    kb = k;
    const IndexType rk = offset + kb;
//
//  Apply the block reflector to the rest of the matrix:
//  A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) -
//                        A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)**T.
//
    if (kb<min(n,m-offset)) {
        blas::mm(NoTrans, Trans,
                 -One, A(_(rk+1,m),_(1,kb)), F(_(kb+1,n),_(1,kb)),
                 One, A(_(rk+1,m),_(kb+1,n)));
    }
//
//  Recomputation of difficult columns.
//
    while (lasticc>0) {
        IndexType iTmp = nint(vn2(lasticc));
        vn1(lasticc) = blas::nrm2(A(_(rk+1,m),lasticc));
//
//      NOTE: The computation of VN1( LSTICC ) relies on the fact that
//      SNRM2 does not fail on vectors with norm below the value of
//      SQRT(DLAMCH('S'))
//
        vn2(lasticc) = vn1(lasticc);
        lasticc = iTmp;
    }
}

//
//  Complex variant
//
template <typename MA, typename JPIV, typename VTAU,
          typename VN1, typename VN2, typename VAUX,
          typename MF>
typename RestrictTo<IsComplex<typename MA::ElementType>::value,
         void>::Type
laqps_impl(typename GeMatrix<MA>::IndexType  offset,
           typename GeMatrix<MA>::IndexType  nb,
           typename GeMatrix<MA>::IndexType  &kb,
           GeMatrix<MA>                      &A,
           DenseVector<JPIV>                 &jPiv,
           DenseVector<VTAU>                 &tau,
           DenseVector<VN1>                  &vn1,
           DenseVector<VN2>                  &vn2,
           DenseVector<VAUX>                 &aux,
           GeMatrix<MF>                      &F)
{
    using std::abs;
    using std::max;
    using std::min;
    using flens::pow;
    using std::sqrt;
    using std::swap;

    typedef typename GeMatrix<MA>::ElementType                  ElementType;
    typedef typename ComplexTrait<ElementType>::PrimitiveType   PrimitiveType;
    typedef typename GeMatrix<MA>::IndexType                    IndexType;

    const Underscore<IndexType> _;

    const IndexType m     = A.numRows();
    const IndexType n     = A.numCols();

    const IndexType lastRk = min(m, n+offset);

    IndexType lasticc = 0;
    IndexType k = 0;

    const PrimitiveType  Zero(0), One(1);
    const ElementType    CZero(0), COne(1);
    const PrimitiveType  tol3z = sqrt(lamch<PrimitiveType>(Eps));
//
//  Beginning of while loop.
//
    while (k<nb && lasticc==0) {
        ++k;

        const IndexType rk = offset + k;
//
//      Determine ith pivot column and swap if necessary
//
        IndexType pvt = (k-1) + blas::iamax(vn1(_(k,n)));
        if (pvt!=k) {
            blas::swap(A(_,pvt), A(_,k));
            blas::swap(F(pvt,_(1,k-1)), F(k,_(1,k-1)));
            swap(jPiv(pvt),jPiv(k));
            vn1(pvt) = vn1(k);
            vn2(pvt) = vn2(k);
        }
//
//      Apply previous Householder reflectors to column K:
//      A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)**T.
//
        if (k>1) {
            blas::conj(F(k,_(1,k-1)));
            blas::mv(NoTrans, -COne, A(_(rk,m),_(1,k-1)), F(k,_(1,k-1)),
                     COne, A(_(rk,m),k));
            blas::conj(F(k,_(1,k-1)));
        }
//
//      Generate elementary reflector H(k).
//
        if (rk<m) {
            larfg(m-rk+1, A(rk,k), A(_(rk+1,m),k), tau(k));
        } else {
            larfg(1, A(m,k), A(_(m+1,m),k), tau(k));
        }

        const ElementType Akk = A(rk,k);
        A(rk,k) = COne;
//
//      Compute Kth column of F:
//
//      Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**T*A(RK:M,K).
//
        if (k<n) {
            blas::mv(ConjTrans, tau(k), A(_(rk,m),_(k+1,n)), A(_(rk,m),k),
                     Zero, F(_(k+1,n),k));
        }
//
//      Padding F(1:K,K) with zeros.
//
        F(_(1,k),k) = CZero;
//
//      Incremental updating of F:
//      F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)**T
//                   *A(RK:M,K).
//
        if (k>1) {
            blas::mv(ConjTrans, -tau(k), A(_(rk,m),_(1,k-1)), A(_(rk,m),k),
                     Zero, aux(_(1,k-1)));
            blas::mv(NoTrans, One, F(_,_(1,k-1)), aux(_(1,k-1)),
                     One, F(_,k));
        }
//
//      Update the current row of A:
//      A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)**H.
//
        if (k<n) {
            blas::mm(NoTrans, ConjTrans,
                     -One, F(_(k+1,n),_(1,k)), A(_(1,rk),_(1,k)),
                     One, A(_(1,rk),_(k+1,n)));
        }
//
//      Update partial column norms.
//
        if (rk<lastRk) {
            for (IndexType j=k+1; j<=n; ++j) {
                if (vn1(j)!=Zero) {
//
//                  NOTE: The following 4 lines follow from the analysis in
//                  Lapack Working Note 176.
//
                    PrimitiveType tmp = abs(A(rk,j)) / vn1(j);
                    tmp = max(Zero, (One+tmp)*(One-tmp));
                    PrimitiveType tmp2 = tmp*pow(vn1(j)/vn2(j),2);
                    if (tmp2<=tol3z) {
                        vn2(j) = PrimitiveType(lasticc);
                        lasticc = j;
                    } else {
                        vn1(j) *= sqrt(tmp);
                    }
                }
            }
        }

        A(rk,k) = Akk;
//
//      End of while loop.
//
    }
    kb = k;
    const IndexType rk = offset + kb;
//
//  Apply the block reflector to the rest of the matrix:
//  A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) -
//                        A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)**T.
//
    if (kb<min(n,m-offset)) {
        blas::mm(NoTrans, ConjTrans,
                 -COne, A(_(rk+1,m),_(1,kb)), F(_(kb+1,n),_(1,kb)),
                 COne, A(_(rk+1,m),_(kb+1,n)));
    }
//
//  Recomputation of difficult columns.
//
    while (lasticc>0) {
        IndexType iTmp = nint(vn2(lasticc));
        vn1(lasticc) = blas::nrm2(A(_(rk+1,m),lasticc));
//
//      NOTE: The computation of VN1( LSTICC ) relies on the fact that
//      SNRM2 does not fail on vectors with norm below the value of
//      SQRT(DLAMCH('S'))
//
        vn2(lasticc) = vn1(lasticc);
        lasticc = iTmp;
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename JPIV, typename VTAU,
          typename VN1, typename VN2, typename VAUX,
          typename MF>
void
laqps_impl(typename GeMatrix<MA>::IndexType  offset,
           typename GeMatrix<MA>::IndexType  nb,
           typename GeMatrix<MA>::IndexType  &kb,
           GeMatrix<MA>                      &A,
           DenseVector<JPIV>                 &jPiv,
           DenseVector<VTAU>                 &tau,
           DenseVector<VN1>                  &vn1,
           DenseVector<VN2>                  &vn2,
           DenseVector<VAUX>                 &aux,
           GeMatrix<MF>                      &F)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    cxxlapack::laqps<IndexType>(A.numRows(),
                                A.numCols(),
                                offset,
                                nb,
                                kb,
                                A.data(),
                                A.leadingDimension(),
                                jPiv.data(),
                                tau.data(),
                                vn1.data(),
                                vn2.data(),
                                aux.data(),
                                F.data(),
                                F.leadingDimension());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
//
//  Real and complex variant
//
template <typename IndexType, typename MA, typename JPIV, typename VTAU,
          typename VN1, typename VN2, typename VAUX, typename MF>
typename RestrictTo<((IsRealGeMatrix<MA>::value &&
                      IsRealDenseVector<VTAU>::value &&
                      IsRealDenseVector<VAUX>::value &&
                      IsRealGeMatrix<MF>::value)
                  || (IsComplexGeMatrix<MA>::value &&
                      IsComplexDenseVector<VTAU>::value &&
                      IsComplexDenseVector<VAUX>::value &&
                      IsComplexGeMatrix<MF>::value))
                  && IsIntegerDenseVector<JPIV>::value
                  && IsRealDenseVector<VN1>::value
                  && IsRealDenseVector<VN2>::value,
         void>::Type
laqps(IndexType   offset,
      IndexType   nb,
      IndexType   &kb,
      MA          &&A,
      JPIV        &&jPiv,
      VTAU        &&tau,
      VN1         &&vn1,
      VN2         &&vn2,
      VAUX        &&aux,
      MF          &&F)
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
    typedef typename RemoveRef<VAUX>::Type      VectorAux;
    typedef typename RemoveRef<MF>::Type        MatrixF;
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
    ASSERT(aux.firstIndex()==1);
    ASSERT(F.firstRow()==1);
    ASSERT(F.firstCol()==1);

    const IndexType n = A.numCols();

    ASSERT(jPiv.length()==n);
    ASSERT(tau.length()>=nb);
    ASSERT(vn1.length()==n);
    ASSERT(vn1.length()==n);
    ASSERT(aux.length()==nb);
    ASSERT(F.numRows()==n);
    ASSERT(F.numCols()==nb);
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
    typename VectorAux::NoView   aux_org    = aux;
    typename MatrixF::NoView     F_org      = F;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::laqps_impl(offset, nb, kb, A, jPiv, tau, vn1, vn2, aux, F);

#   ifdef CHECK_CXXLAPACK
//
//  Restore output arguments
//
    typename MatrixA::NoView     A_generic    = A;
    typename VectorJPiv::NoView  jPiv_generic = jPiv;
    typename VectorTau::NoView   tau_generic  = tau;
    typename VectorVN1::NoView   vn1_generic  = vn1;
    typename VectorVN2::NoView   vn2_generic  = vn2;
    typename VectorAux::NoView   aux_generic  = aux;
    typename MatrixF::NoView     F_generic    = F;

    A    = A_org;
    jPiv = jPiv_org;
    tau  = tau_org;
    vn1  = vn1_org;
    vn2  = vn2_org;
    aux  = aux_org;
    F    = F_org;

//
//  Compare results
//
    external::laqps_impl(offset, nb, kb, A, jPiv, tau, vn1, vn2, aux, F);

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

    if (! isIdentical(aux_generic, aux, "aux_generic", "aux")) {
        std::cerr << "CXXLAPACK: aux_generic = " << aux_generic << std::endl;
        std::cerr << "F77LAPACK: aux = " << aux << std::endl;
        failed = true;
    }

    if (! isIdentical(F_generic, F, "F_generic", "F")) {
        std::cerr << "CXXLAPACK: F_generic = " << F_generic << std::endl;
        std::cerr << "F77LAPACK: F = " << F << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAQPS_TCC
