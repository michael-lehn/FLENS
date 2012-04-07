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
       SUBROUTINE DGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1)                                  --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2009                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_HRD_TCC
#define FLENS_LAPACK_IMPL_HRD_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename IndexType, typename MA>
IndexType
hrd_wsq_impl(IndexType iLo, IndexType iHi, const GeMatrix<MA> &A)
{
    using std::min;

    typedef typename GeMatrix<MA>::ElementType  T;

//
//  Perform and apply workspace query
//
    const IndexType nbMax = 64;
    const IndexType n   = A.numCols();

    return n*min(nbMax, IndexType(ilaenv<T>(1,"GEHRD", "", n, iLo, iHi)));
}

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
hrd_impl(IndexType           iLo,
         IndexType           iHi,
         GeMatrix<MA>        &A,
         DenseVector<VTAU>   &tau,
         DenseVector<VWORK>  &work)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::View         GeView;
    typedef typename GeView::Engine             GeViewEngine;

    const Underscore<IndexType> _;

//
//  Paramters, e.g. maximum block size and buffer on stack for TrMatrix Tr.
//
    const T             One(1);
    const IndexType     nbMax = 64;
    const IndexType     ldt = nbMax + 1;
    T                   trBuffer[nbMax*ldt];
    GeView              _Tr = GeViewEngine(ldt, nbMax, trBuffer, ldt);
//
//  Perform and apply workspace query
//
    const IndexType n = A.numCols();
    IndexType       nb = min(nbMax,
                             IndexType(ilaenv<T>(1,"GEHRD", "", n, iLo, iHi)));

    const IndexType lWorkOpt = n*nb;
    if (work.length()==0) {
        work.resize(lWorkOpt);
        work(1) = lWorkOpt;
    }
//
//  Set elements 1:iLo-1 and iHi:n-1 of tau to zero
//
    tau(_(1,iLo-1)) = 0;
    tau(_(max(IndexType(1),iHi),n-1)) = 0;
//
//  Quick return if possible
//
    const IndexType nh = iHi - iLo + 1;
    if (nh<=1) {
        work(1) = 1;
        return;
    }
//
//  Determine the block size
//
    IndexType nbMin = 2;
    IndexType iws = 1;
    IndexType nx = 0;

    if (nb>1 && nb<nh) {
//
//      Determine when to cross over from blocked to unblocked code
//      (last block is always handled by unblocked code)
//
        nx = max(nb, IndexType(ilaenv<T>(3, "GEHRD", "", n, iLo, iHi)));
        if (nx<nh) {
//
//          Determine if workspace is large enough for blocked code
//
            iws = n*nb;
            if (work.length()<iws) {
//
//              Not enough workspace to use optimal NB:  determine the
//              minimum value of NB, and reduce NB or force use of
//              unblocked code
//
                nbMin = max(2, ilaenv<T>(2, "GEHRD", "", n, iLo, iHi));
                if (work.length()>=n*nbMin) {
                    nb = work.length() / n;
                } else {
                    nb = 1;
                }
            }
        }
    }
    IndexType ldWork = n;
    IndexType i = iLo;
    if (nb<nbMin || nb>=nh) {
//
//      Use unblocked code below
//
    } else {
//
//      Use blocked code
//
        for (i=iLo; i<=iHi-1-nx; i+=nb) {
            const IndexType ib = min(nb, iHi-i);
            auto  Tr = _Tr(_(1,ib),_(1,ib)).upper();
            GeView  Y(iHi, ib, work, ldWork);
//
//          Reduce columns i:i+ib-1 to Hessenberg form, returning the
//          matrices V and T of the block reflector H = I - V*T*V**T
//          which performs the reduction, and also the matrix Y = A*V*T
//
            lahr2(i, ib, A(_(1,iHi),_(i,iHi)), tau(_(i,i+ib-1)), Tr, Y);
//
//          Partition V and the parts of matrix A that get modified
//
            auto V1 = A(_(i+1, i+ib-1),_(i,i+ib-2));
            auto V2 = A(_(i+ib,   iHi),_(i,i+ib-1));
//
//          Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
//          right, computing  A := A - Y * V**T. V(i+ib,ib-1) must be set
//          to 1
//
            T ei = A(i+ib,i+ib-1);
            A(i+ib,i+ib-1) = One;
            blas::mm(NoTrans, Trans, -One, Y, V2, One, A(_(1,iHi),_(i+ib,iHi)));
            A(i+ib,i+ib-1) = ei;
//
//          Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
//          right
//
            blas::mm(Right, Trans, One, V1.lowerUnit(), Y(_(1,i),_(1,ib-1)));
            A(_(1,i),_(i+1,i+ib-1)) -= Y(_(1,i),_(1,ib-1));
//
//          Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
//          left
//
            GeView  Work(n-i-ib+1, ib, work, ldWork);
            larfb(Left, Trans, Forward, ColumnWise,
                  A(_(i+1,iHi),_(i,i+ib-1)),
                  Tr,
                  A(_(i+1,iHi),_(i+ib,n)),
                  Work);
        }
    }
//
//  Use unblocked code to reduce the rest of the matrix
//
    hd2(i, iHi, A, tau, work);
    work(1) = iws;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename IndexType, typename MA>
IndexType
hrd_wsq_impl(IndexType            iLo,
             IndexType            iHi,
             const GeMatrix<MA>   &A)
{
    typedef typename GeMatrix<MA>::ElementType  T;

    T                   WORK, DUMMY;
    const IndexType     LWORK  = -1;

    cxxlapack::gehrd<IndexType>(A.numRows(),
                                iLo,
                                iHi,
                                &DUMMY,
                                A.leadingDimension(),
                                &DUMMY,
                                &WORK,
                                LWORK);
    return WORK;
}

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
hrd_impl(IndexType            iLo,
         IndexType            iHi,
         GeMatrix<MA>         &A,
         DenseVector<VTAU>    &tau,
         DenseVector<VWORK>   &work)
{
    cxxlapack::gehrd<IndexType>(A.numRows(),
                                iLo,
                                iHi,
                                A.data(),
                                A.leadingDimension(),
                                tau.data(),
                                work.data(),
                                work.length());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename IndexType, typename MA>
IndexType
hrd_wsq(IndexType           iLo,
        IndexType           iHi,
        const GeMatrix<MA>  &A)
{
    LAPACK_DEBUG_OUT("hrd_wsq");

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::hrd_wsq_impl(iLo, iHi, A);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    IndexType _info = external::hrd_wsq_impl(iLo, iHi, A);

    if (! isIdentical(info, _info, " info", "_info")) {
        ASSERT(0);
    }
#   endif

    return info;
}

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
hrd(IndexType           iLo,
    IndexType           iHi,
    GeMatrix<MA>        &A,
    DenseVector<VTAU>   &tau,
    DenseVector<VWORK>  &work)
{
    LAPACK_DEBUG_OUT("hrd");

    using std::max;
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());
    ASSERT(tau.firstIndex()==1);

    const IndexType n = A.numCols();

    if (n==0) {
        ASSERT(iLo==1);
        ASSERT(iHi==0);
    } else {
        ASSERT(1<=iLo);
        ASSERT(iLo<=iHi);
        ASSERT(iHi<=n);
    }

    ASSERT(tau.length()==(n-1));
    const bool lQuery = (work.length()==0);
    ASSERT(lQuery || (work.length()>=max(IndexType(1),n)));
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename GeMatrix<MA>::NoView       _A = A;
    typename DenseVector<VTAU>::NoView  _tau = tau;
    typename DenseVector<VWORK>::NoView _work = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::hrd_impl(iLo, iHi, A, tau, work);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    if (_work.length()==0) {
        _work.resize(work.length());
    }
    external::hrd_impl(iLo, iHi, _A, _tau, _work);

    bool failed = false;
    if (! isIdentical(A, _A, " A", "_A")) {
        std::cerr << "CXXLAPACK:  A = " << A << std::endl;
        std::cerr << "F77LAPACK: _A = " << _A << std::endl;
        failed = true;
    }

    if (! isIdentical(tau, _tau, " tau", "_tau")) {
        std::cerr << "CXXLAPACK:  tau = " << tau << std::endl;
        std::cerr << "F77LAPACK: _tau = " << _tau << std::endl;
        failed = true;
    }

    if (! isIdentical(work, _work, " work", "_work")) {
        std::cerr << "CXXLAPACK:  work = " << work << std::endl;
        std::cerr << "F77LAPACK: _work = " << _work << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: hrd.tcc" << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: hrd.tcc" << std::endl;
    }

#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename MA>
IndexType
hrd_wsq(IndexType   iLo,
        IndexType   iHi,
        const MA    &&A)
{
    CHECKPOINT_ENTER;
    const IndexType info = hrd_wsq(iLo, iHi, A);
    CHECKPOINT_LEAVE;

    return info;
}

template <typename IndexType, typename MA,
          typename VTAU, typename VWORK>
void
hrd(IndexType       iLo,
    IndexType       iHi,
    MA              &&A,
    VTAU            &&tau,
    VWORK           &&work)
{
    CHECKPOINT_ENTER;
    hrd(iLo, iHi, A, tau, work);
    CHECKPOINT_LEAVE;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_HRD_TCC
