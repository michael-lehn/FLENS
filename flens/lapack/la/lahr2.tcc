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
       SUBROUTINE DLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
       SUBROUTINE ZLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
 *
 *  -- LAPACK auxiliary routine (version 3.3.1)                        --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2009                                                      --
 */

#ifndef FLENS_LAPACK_LA_LAHR2_TCC
#define FLENS_LAPACK_LA_LAHR2_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//
//  Real variant
//
template <typename IndexType, typename MA, typename VTAU,
          typename MTR, typename MY>
typename RestrictTo<IsReal<typename MA::ElementType>::value,
         void>::Type
lahr2_impl(IndexType k, IndexType nb, GeMatrix<MA> &A,
           DenseVector<VTAU> &tau,   TrMatrix<MTR> &Tr,
           GeMatrix<MY> &Y)
{
    using lapack::larfg;
    using std::min;

    typedef typename GeMatrix<MA>::ElementType  T;

    const Underscore<IndexType> _;

    const IndexType n = A.numRows();
    const T         Zero(0), One(1);

//  TODO: as long as view creation is not supported by the TrMatrix interface
//        we get them through a GeMatrix::View
    auto _Tr = Tr.general();
//
//  Quick return if possible
//
    if (n<=1) {
        return;
    }

    T ei = T(0);
    for (IndexType i=1; i<=nb; ++i) {
        if (i>1) {
//
//          Update A(K+1:N,I)
//
//          Update I-th column of A - Y * V**T
//
            blas::mv(NoTrans,
                     -One, Y(_(k+1,n),_(1,i-1)), A(k+i-1,_(1,i-1)),
                     One, A(_(k+1,n),i));
//
//          Apply I - V * T**T * V**T to this column (call it b) from the
//          left, using the last column of T as workspace
//
//          Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
//                   ( V2 )             ( b2 )
//
//          where V1 is unit lower triangular
//
//          w := V1**T * b1
//
            _Tr(_(1,i-1),nb) = A(_(k+1,k+i-1),i);
            blas::mv(Trans, A(_(k+1,k+i-1),_(1,i-1)).lowerUnit(),
                     _Tr(_(1,i-1),nb));
//
//          w := w + V2**T * b2
//
            blas::mv(Trans,
                     One, A(_(k+i,n),_(1,i-1)), A(_(k+i,n),i),
                     One, _Tr(_(1,i-1),nb));
//
//          w := T**T * w
//
            blas::mv(Trans, _Tr(_(1,i-1),_(1,i-1)).upper(), _Tr(_(1,i-1),nb));
//
//          b2 := b2 - V2*w
//
            blas::mv(NoTrans,
                     -One, A(_(k+i,n),_(1,i-1)), _Tr(_(1,i-1),nb),
                     One, A(_(k+i,n),i));
//
//          b1 := b1 - V1*w
//
            blas::mv(NoTrans,
                     A(_(k+1,k+i-1),_(1,i-1)).lowerUnit(),
                     _Tr(_(1,i-1),nb));
            A(_(k+1,k+i-1),i) -= _Tr(_(1,i-1),nb);

            A(k+i-1,i-1) = ei;
        }
//
//      Generate the elementary reflector H(I) to annihilate
//      A(K+I+1:N,I)
//
        larfg(n-k-i+1, A(k+i,i), A(_(min(k+i+1,n),n), i), tau(i));

        ei = A(k+i, i);
        A(k+i, i) = One;
//
//      Compute  Y(K+1:N,I)
//
        blas::mv(NoTrans,
                 One, A(_(k+1,n),_(i+1,n-k+1)), A(_(k+i,n),i),
                 Zero, Y(_(k+1,n),i));
        blas::mv(Trans,
                 One, A(_(k+i,n),_(1,i-1)), A(_(k+i,n),i),
                 Zero, _Tr(_(1,i-1),i));
        blas::mv(NoTrans,
                 -One, Y(_(k+1,n),_(1,i-1)), _Tr(_(1,i-1),i),
                 One, Y(_(k+1,n),i));
        blas::scal(tau(i),Y(_(k+1,n),i));
//
//      Compute T(1:I,I)
//
        blas::scal(-tau(i), _Tr(_(1,i-1),i));
        blas::mv(NoTrans, _Tr(_(1,i-1),_(1,i-1)).upper(), _Tr(_(1,i-1),i));
        Tr(i,i) = tau(i);
    }
    A(k+nb, nb) = ei;
//
//  Compute Y(1:K,1:NB)
//
    Y(_(1,k),_(1,nb)) = A(_(1,k),_(2,2+nb-1));
    blas::mm(Right, NoTrans,
             One, A(_(k+1,k+nb),_(1,nb)).lowerUnit(),
             Y(_(1,k),_(1,nb)));

    if (n>k+nb) {
        blas::mm(NoTrans, NoTrans,
                 One, A(_(1,k),_(2+nb,n-k+1)), A(_(k+1+nb,n),_(1,nb)),
                 One, Y(_(1,k),_(1,nb)));
    }
    blas::mm(Right, NoTrans, One, Tr, Y(_(1,k),_(1,nb)));
}

//
//  Complex variant
//
template <typename IndexType, typename MA, typename VTAU,
          typename MTR, typename MY>
typename RestrictTo<IsComplex<typename MA::ElementType>::value,
         void>::Type
lahr2_impl(IndexType k, IndexType nb, GeMatrix<MA> &A,
           DenseVector<VTAU> &tau,   TrMatrix<MTR> &Tr,
           GeMatrix<MY> &Y)
{
    using lapack::larfg;
    using std::min;

    typedef typename GeMatrix<MA>::ElementType  T;

    const Underscore<IndexType> _;

    const IndexType n = A.numRows();
    const T         Zero(0), One(1);

//  TODO: as long as view creation is not supported by the TrMatrix interface
//        we get them through a GeMatrix::View
    auto _Tr = Tr.general();
//
//  Quick return if possible
//
    if (n<=1) {
        return;
    }

    T ei = Zero;
    for (IndexType i=1; i<=nb; ++i) {
        if (i>1) {
//
//          Update A(K+1:N,I)
//
//          Update I-th column of A - Y * V**H
//
            blas::conj(A(k+i-1,_(1,i-1)));
            blas::mv(NoTrans,
                     -One, Y(_(k+1,n),_(1,i-1)), A(k+i-1,_(1,i-1)),
                     One, A(_(k+1,n),i));
            blas::conj(A(k+i-1,_(1,i-1)));
//
//          Apply I - V * T**H * V**H to this column (call it b) from the
//          left, using the last column of T as workspace
//
//          Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
//                   ( V2 )             ( b2 )
//
//          where V1 is unit lower triangular
//
//          w := V1**H * b1
//
            _Tr(_(1,i-1),nb) = A(_(k+1,k+i-1),i);
            blas::mv(ConjTrans, A(_(k+1,k+i-1),_(1,i-1)).lowerUnit(),
                     _Tr(_(1,i-1),nb));
//
//          w := w + V2**H * b2
//
            blas::mv(ConjTrans,
                     One, A(_(k+i,n),_(1,i-1)), A(_(k+i,n),i),
                     One, _Tr(_(1,i-1),nb));
//
//          w := T**H * w
//
            blas::mv(ConjTrans, _Tr(_(1,i-1),_(1,i-1)).upper(),
                     _Tr(_(1,i-1),nb));
//
//          b2 := b2 - V2*w
//
            blas::mv(NoTrans,
                     -One, A(_(k+i,n),_(1,i-1)), _Tr(_(1,i-1),nb),
                     One, A(_(k+i,n),i));
//
//          b1 := b1 - V1*w
//
            blas::mv(NoTrans,
                     A(_(k+1,k+i-1),_(1,i-1)).lowerUnit(),
                     _Tr(_(1,i-1),nb));
            A(_(k+1,k+i-1),i) -= _Tr(_(1,i-1),nb);

            A(k+i-1,i-1) = ei;
        }
//
//      Generate the elementary reflector H(I) to annihilate
//      A(K+I+1:N,I)
//
        larfg(n-k-i+1, A(k+i,i), A(_(min(k+i+1,n),n), i), tau(i));

        ei = A(k+i, i);
        A(k+i, i) = One;
//
//      Compute  Y(K+1:N,I)
//
        blas::mv(NoTrans,
                 One, A(_(k+1,n),_(i+1,n-k+1)), A(_(k+i,n),i),
                 Zero, Y(_(k+1,n),i));
        blas::mv(ConjTrans,
                 One, A(_(k+i,n),_(1,i-1)), A(_(k+i,n),i),
                 Zero, _Tr(_(1,i-1),i));
        blas::mv(NoTrans,
                 -One, Y(_(k+1,n),_(1,i-1)), _Tr(_(1,i-1),i),
                 One, Y(_(k+1,n),i));
        blas::scal(tau(i),Y(_(k+1,n),i));
//
//      Compute T(1:I,I)
//
        blas::scal(-tau(i), _Tr(_(1,i-1),i));
        blas::mv(NoTrans, _Tr(_(1,i-1),_(1,i-1)).upper(), _Tr(_(1,i-1),i));
        Tr(i,i) = tau(i);
    }
    A(k+nb, nb) = ei;
//
//  Compute Y(1:K,1:NB)
//
    Y(_(1,k),_(1,nb)) = A(_(1,k),_(2,2+nb-1));
    blas::mm(Right, NoTrans,
             One, A(_(k+1,k+nb),_(1,nb)).lowerUnit(),
             Y(_(1,k),_(1,nb)));

    if (n>k+nb) {
        blas::mm(NoTrans, NoTrans,
                 One, A(_(1,k),_(2+nb,n-k+1)), A(_(k+1+nb,n),_(1,nb)),
                 One, Y(_(1,k),_(1,nb)));
    }
    blas::mm(Right, NoTrans, One, Tr, Y(_(1,k),_(1,nb)));
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename IndexType, typename MA, typename VTAU,
          typename MTR, typename MY>
void
lahr2_impl(IndexType k, IndexType nb, GeMatrix<MA> &A,
           DenseVector<VTAU> &tau,   TrMatrix<MTR> &Tr,
           GeMatrix<MY> &Y)
{
    cxxlapack::lahr2<IndexType>(A.numRows(),
                                k,
                                nb,
                                A.data(),
                                A.leadingDimension(),
                                tau.data(),
                                Tr.data(),
                                Tr.leadingDimension(),
                                Y.data(),
                                Y.leadingDimension());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
//
//  Real variant
//
template <typename IndexType, typename MA, typename VTAU,
          typename MTR, typename MY>
typename RestrictTo<IsRealGeMatrix<MA>::value
                 && IsRealDenseVector<VTAU>::value
                 && IsRealTrMatrix<MTR>::value
                 && IsRealGeMatrix<MY>::value,
         void>::Type
lahr2(IndexType     k,
      IndexType     nb,
      MA            &&A,
      VTAU          &&tau,
      MTR           &&Tr,
      MY            &&Y)
{
    LAPACK_DEBUG_OUT("lahr2");

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MA>::Type        MatrixA;
    typedef typename RemoveRef<VTAU>::Type      VectorTau;
    typedef typename RemoveRef<MTR>::Type       MatrixTr;
    typedef typename RemoveRef<MY>::Type        MatrixY;
#   endif

//
//  Test the input parameters
//
    ASSERT(k<A.numRows());
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numCols()==A.numRows()-k+1);
    ASSERT(tau.length()==nb);
    ASSERT(Tr.dim()==nb);
    ASSERT(Y.numRows()==A.numRows());
    ASSERT(Y.numCols()==nb);

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView           A_org     = A;
    typename VectorTau::NoView         tau_org   = tau;
    typename MatrixTr::GeneralNoView   Tr_org    = Tr.general();
    typename MatrixY::NoView           Y_org     = Y;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::lahr2_impl(k, nb, A, tau, Tr, Y);

#   ifdef CHECK_CXXLAPACK
//
//  Restore output arguments
//
    typename MatrixA::NoView           A_generic     = A;
    typename VectorTau::NoView         tau_generic   = tau;
    typename MatrixTr::GeneralNoView   Tr_generic    = Tr.general();
    typename MatrixY::NoView           Y_generic     = Y;

    A            = A_org;
    tau          = tau_org;
    Tr.general() = Tr_org;
    Y            = Y_org;

//
//  Compare results
//
    external::lahr2_impl(k, nb, A, tau, Tr, Y);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(tau_generic, tau, "tau_generic", "tau")) {
        std::cerr << "CXXLAPACK: tau_generic = " << tau_generic << std::endl;
        std::cerr << "F77LAPACK: tau = " << tau << std::endl;
        failed = true;
    }

    if (! isIdentical(Tr_generic, Tr.general(), " Tr_generic", "Tr")) {
        std::cerr << "CXXLAPACK: Tr_generic = " << Tr_generic << std::endl;
        std::cerr << "F77LAPACK: Tr = " << Tr.general() << std::endl;
        failed = true;
    }

    if (! isIdentical(Y_generic, Y, "Y_generic", "Y")) {
        std::cerr << "CXXLAPACK: Y_generic = " << Y_generic << std::endl;
        std::cerr << "F77LAPACK: Y = " << Y << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: lahr2.tcc" << std::endl;
        std::cerr << "k = " << k << ", nb = " << nb << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: lahr2.tcc" << std::endl;
    }
#   endif
}

//
//  Complex variant
//
template <typename IndexType, typename MA, typename VTAU,
          typename MTR, typename MY>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexDenseVector<VTAU>::value
                 && IsComplexTrMatrix<MTR>::value
                 && IsComplexGeMatrix<MY>::value,
         void>::Type
lahr2(IndexType     k,
      IndexType     nb,
      MA            &&A,
      VTAU          &&tau,
      MTR           &&Tr,
      MY            &&Y)
{
    LAPACK_DEBUG_OUT("lahr2");

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MA>::Type        MatrixA;
    typedef typename RemoveRef<VTAU>::Type      VectorTau;
    typedef typename RemoveRef<MTR>::Type       MatrixTr;
    typedef typename RemoveRef<MY>::Type        MatrixY;
#   endif

//
//  Test the input parameters
//
    ASSERT(k<A.numRows());
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numCols()==A.numRows()-k+1);
    ASSERT(tau.length()==nb);
    ASSERT(Tr.dim()==nb);
    ASSERT(Y.numRows()==A.numRows());
    ASSERT(Y.numCols()==nb);

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView           A_org     = A;
    typename VectorTau::NoView         tau_org   = tau;
    typename MatrixTr::GeneralNoView   Tr_org    = Tr.general();
    typename MatrixY::NoView           Y_org     = Y;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::lahr2_impl(k, nb, A, tau, Tr, Y);

#   ifdef CHECK_CXXLAPACK
//
//  Restore output arguments
//
    typename MatrixA::NoView           A_generic     = A;
    typename VectorTau::NoView         tau_generic   = tau;
    typename MatrixTr::GeneralNoView   Tr_generic    = Tr.general();
    typename MatrixY::NoView           Y_generic     = Y;

    A            = A_org;
    tau          = tau_org;
    Tr.general() = Tr_org;
    Y            = Y_org;

//
//  Compare results
//
    external::lahr2_impl(k, nb, A, tau, Tr, Y);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(tau_generic, tau, "tau_generic", "tau")) {
        std::cerr << "CXXLAPACK: tau_generic = " << tau_generic << std::endl;
        std::cerr << "F77LAPACK: tau = " << tau << std::endl;
        failed = true;
    }

    if (! isIdentical(Tr_generic, Tr.general(), " Tr_generic", "Tr")) {
        std::cerr << "CXXLAPACK: Tr_generic = " << Tr_generic << std::endl;
        std::cerr << "F77LAPACK: Tr = " << Tr.general() << std::endl;
        failed = true;
    }

    if (! isIdentical(Y_generic, Y, "Y_generic", "Y")) {
        std::cerr << "CXXLAPACK: Y_generic = " << Y_generic << std::endl;
        std::cerr << "F77LAPACK: Y = " << Y << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: lahr2.tcc" << std::endl;
        std::cerr << "k = " << k << ", nb = " << nb << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: lahr2.tcc" << std::endl;
    }
#   endif
}


} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAHR2_TCC
