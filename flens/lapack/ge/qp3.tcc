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
       SUBROUTINE DGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO )
       SUBROUTINE ZGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, RWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_GE_QP3_TCC
#define FLENS_LAPACK_GE_QP3_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- worksize query ------------------------------------------------------------

template <typename MA>
typename GeMatrix<MA>::IndexType
qp3_wsq_impl(const GeMatrix<MA> &A)
{
    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const IndexType m     = A.numRows();
    const IndexType n     = A.numCols();
    const IndexType minmn = min(m, n);

    IndexType iws, lwOpt;

    if (minmn==0) {
        iws   = 1;
        lwOpt = 1;
    } else {
        iws = 3*n + 1;
        const IndexType nb = ilaenv<ElementType>(1, "GEQRF", "", m, n);
        lwOpt = 2*n +(n+1)*nb;
    }

    return lwOpt;
}


//-- (ge)qp3 [real variant] ----------------------------------------------------

template <typename MA, typename JPIV, typename VTAU, typename VWORK>
void
qp3_impl(GeMatrix<MA> &A, DenseVector<JPIV> &jPiv, DenseVector<VTAU> &tau,
         DenseVector<VWORK> &work)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const Underscore<IndexType> _;

    const IndexType m     = A.numRows();
    const IndexType n     = A.numCols();
    const IndexType minmn = min(m, n);

    IndexType iws, lwOpt;

    if (minmn==0) {
        iws   = 1;
    } else {
        iws = 3*n + 1;
    }

    lwOpt = qp3_wsq(A);

    if (work.length()==0) {
        work.resize(lwOpt);
    }
    work(1) = lwOpt;

    IndexType lWork = work.length();

//
//  Quick return if possible.
//
    if (minmn==0) {
        return;
    }
//
//  Move initial columns up front.
//
    IndexType  nFixed = 1;
    for (IndexType j=1; j<=n; ++j) {
        if (jPiv(j)!=0) {
            if (j!=nFixed) {
                blas::swap(A(_,j), A(_,nFixed));
                jPiv(j)      = jPiv(nFixed);
                jPiv(nFixed) = j;
            } else {
                jPiv(j) = j;
            }
            ++nFixed;
        } else {
            jPiv(j) = j;
        }
    }
    --nFixed;
//
//  Factorize fixed columns
//  =======================
//
//  Compute the QR factorization of fixed columns and update
//  remaining columns.
//
    if (nFixed>0) {
        const IndexType na = min(m, nFixed);
        auto A1 = A(_,_(1,na));
        auto A2 = A(_,_(na+1,n));
        auto tau1 = tau(_(1,na));

        qrf(A1, tau1, work);
        iws = max(iws, IndexType(work(1)));
        if (na<n) {
            ormqr(Left, Trans, A1, tau1, A2, work);
            iws = max(iws, IndexType(work(1)));
        }
    }
//
//  Factorize free columns
//  ======================
//
    if (nFixed<minmn) {

        IndexType sm = m - nFixed;
        IndexType sn = n - nFixed;
        IndexType sminmn = minmn - nFixed;
//
//      Determine the block size.
//
        IndexType nb = ilaenv<T>(1, "GEQRF", "", sm, sn);
        IndexType nbMin = 2;
        IndexType nx = 0;

        if (nb>1 && nb<sminmn) {
//
//          Determine when to cross over from blocked to unblocked code.
//
            nx = max(IndexType(0), ilaenv<T>(3, "GEQRF", "", sm, sn));
//
//
            if (nx<sminmn) {
//
//              Determine if workspace is large enough for blocked code.
//
                IndexType minWs = 2*sn + (sn+1)*nb;
                iws = max(iws, minWs);
                if (lWork<minWs) {
//
//                  Not enough workspace to use optimal NB: Reduce NB and
//                  determine the minimum value of NB.
//
                    nb = (lWork-2*sn) / (sn+1);
                    nbMin = max(IndexType(2),
                                ilaenv<T>(2, "GEQRF", "", sm, sn));

                }
            }
        }
//
//      Initialize partial column norms. The first N elements of work
//      store the exact column norms.
//
        for (IndexType j=nFixed+1; j<=n; ++j) {
            work(j) = blas::nrm2(A(_(nFixed+1,m),j));
            work(n+j) = work(j);
        }

        IndexType j;

        if (nb>=nbMin && nb<sminmn && nx<sminmn) {
//
//          Use blocked code initially.
//
            j = nFixed + 1;
//
//          Compute factorization: while loop.
//
//
            const IndexType topbmn = minmn - nx;
            while (j<=topbmn) {
                const IndexType jb = min(nb, topbmn-j+1);
//
//              Factorize JB columns among columns J:N.
//
                IndexType fjb;
                auto _A    = A(_,_(j,n));
                auto _jPiv = jPiv(_(j,n));
                auto _tau  = tau(_(j,min(j+jb-1,minmn)));
                auto vn1   = work(_(j,n));
                auto vn2   = work(_(j+n,2*n));
                auto aux   = work(_(2*n+1, 2*n+jb));

                IndexType fLen = jb*(n-j+1);
                auto      f    = work(_(2*n+jb+1,2*n+jb+fLen));

                GeMatrixView<T> F(n-j+1, jb, f, n-j+1);

                laqps(j-1, jb, fjb, _A, _jPiv, _tau, vn1, vn2, aux, F);

                j += fjb;
            }
        } else {
            j = nFixed + 1;
        }
//
//      Use unblocked code to factor the last or only block.
//
//
        if (j<=minmn) {
            auto _A    = A(_,_(j,n));
            auto _jPiv = jPiv(_(j,n));
            auto _tau  = tau(_(j,minmn));
            auto vn1   = work(_(j,n));
            auto vn2   = work(_(j+n,2*n));
            auto _work = work(_(2*n+1, 3*n+1-j));

            laqp2(j-1, _A, _jPiv, _tau, vn1, vn2, _work);
        }

    }

    work(1) = iws;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- worksize query ------------------------------------------------------------

template <typename MA>
typename RestrictTo<IsRealGeMatrix<MA>::value,
         typename MA::IndexType>::Type
qp3_wsq_impl(const MA &A)
{
    typedef typename MA::ElementType  ElementType;
    typedef typename MA::IndexType    IndexType;

    ElementType  DUMMY, WORK;
    IndexType    IDUMMY;
    IndexType    LWORK = -1;

    cxxlapack::geqp3<IndexType>(A.numRows(),
                                A.numCols(),
                                &DUMMY,
                                A.leadingDimension(),
                                &IDUMMY,
                                &DUMMY,
                                &WORK,
                                LWORK);
    return WORK;
}

template <typename MA>
typename RestrictTo<IsComplexGeMatrix<MA>::value,
         typename MA::IndexType>::Type
qp3_wsq_impl(const MA &A)
{
    typedef typename MA::ElementType                            ElementType;
    typedef typename ComplexTrait<ElementType>::PrimitiveType   RealType;
    typedef typename MA::IndexType                              IndexType;

    ElementType  DUMMY, WORK;
    RealType     RWORK;
    IndexType    IDUMMY;
    IndexType    LWORK = -1;

    cxxlapack::geqp3<IndexType>(A.numRows(),
                                A.numCols(),
                                &DUMMY,
                                A.leadingDimension(),
                                &IDUMMY,
                                &DUMMY,
                                &WORK,
                                LWORK,
                                &RWORK);
    return WORK.real();
}

//-- (ge)qp3 [real variant] ----------------------------------------------------

template <typename MA, typename JPIV, typename VTAU, typename VWORK>
void
qp3_impl(GeMatrix<MA>       &A,
         DenseVector<JPIV>  &jPiv,
         DenseVector<VTAU>  &tau,
         DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    cxxlapack::geqp3<IndexType>(A.numRows(),
                                A.numCols(),
                                A.data(),
                                A.leadingDimension(),
                                jPiv.data(),
                                tau.data(),
                                work.data(),
                                work.length());
}

//-- (ge)qp3 [complex variant] -------------------------------------------------

template <typename MA, typename JPIV, typename VTAU, typename VWORK,
          typename VRWORK>
void
qp3_impl(GeMatrix<MA>        &A,
         DenseVector<JPIV>   &jPiv,
         DenseVector<VTAU>   &tau,
         DenseVector<VWORK>  &work,
         DenseVector<VRWORK> &rWork)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    cxxlapack::geqp3<IndexType>(A.numRows(),
                                A.numCols(),
                                A.data(),
                                A.leadingDimension(),
                                jPiv.data(),
                                tau.data(),
                                work.data(),
                                work.length(),
                                rWork.data());
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- (ge)qp3 [real variant] ----------------------------------------------------

template <typename MA, typename VJPIV, typename VTAU, typename VWORK>
typename RestrictTo<IsRealGeMatrix<MA>::value
                 && IsIntegerDenseVector<VJPIV>::value
                 && IsRealDenseVector<VTAU>::value
                 && IsRealDenseVector<VWORK>::value,
         void>::Type
qp3(MA      &&A,
    VJPIV   &&jPiv,
    VTAU    &&tau,
    VWORK   &&work)
{
    using std::min;
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type        MatrixA;
    typedef typename MatrixA::IndexType         IndexType;
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<VJPIV>::Type     VectorJPiv;
    typedef typename RemoveRef<VTAU>::Type      VectorTau;
    typedef typename RemoveRef<VWORK>::Type     VectorWork;
#   endif

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();
    const IndexType k = min(m, n);

#   ifndef NDEBUG
//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(jPiv.firstIndex()==1);
    ASSERT(tau.firstIndex()==1);
    ASSERT(work.firstIndex()==1);

    ASSERT(jPiv.length()==0 || jPiv.length()==n);
    ASSERT(tau.length()==0 || tau.length()==k);
    ASSERT(work.length()>=3*n+1 || work.length()==IndexType(0));
#   endif

//
//  Resize vector output arguments and workspace if needed
//
    if (jPiv.length()==0) {
        jPiv.resize(n);
    }
    if (tau.length()==0) {
        tau.resize(k);
    }
    if (work.length()==0) {
        work.resize(qp3_wsq(A));
    }

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename MatrixA::NoView     A_org      = A;
    typename VectorJPiv::NoView  jPiv_org   = jPiv;
    typename VectorTau::NoView   tau_org    = tau;
    typename VectorWork::NoView  work_org   = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::qp3_impl(A, jPiv, tau, work);

#   ifdef CHECK_CXXLAPACK
//
//  Restore output arguments
//
    typename MatrixA::NoView     A_generic    = A;
    typename VectorJPiv::NoView  jPiv_generic = jPiv;
    typename VectorTau::NoView   tau_generic  = tau;
    typename VectorWork::NoView  work_generic = work;

    A    = A_org;
    jPiv = jPiv_org;
    tau  = tau_org;

    // if the generic implementation resized work due to a work size query
    // we must not restore the work array
    if (work_org.length()>0) {
        work = work_org;
    } else {
        work = 0;
    }
//
//  Compare results
//
    external::qp3_impl(A, jPiv, tau, work);

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

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

//
//-- (ge)qp3 [real variant (with temporary workspace) ] ------------------------
//
template <typename MA, typename VJPIV, typename VTAU>
typename RestrictTo<IsRealGeMatrix<MA>::value
                 && IsIntegerDenseVector<VJPIV>::value
                 && IsRealDenseVector<VTAU>::value,
         void>::Type
qp3(MA      &&A,
    VJPIV   &&jPiv,
    VTAU    &&tau)
{
    typename RemoveRef<MA>::Type::Vector  work;

    qp3(A, jPiv, tau, work);
}


#ifdef USE_CXXLAPACK
//
//-- (ge)qp3 [complex variant] -------------------------------------------------
//
template <typename MA, typename VJPIV, typename VTAU, typename VWORK,
          typename VRWORK>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsIntegerDenseVector<VJPIV>::value
                 && IsComplexDenseVector<VTAU>::value
                 && IsComplexDenseVector<VWORK>::value
                 && IsRealDenseVector<VRWORK>::value,
         void>::Type
qp3(MA      &&A,
    VJPIV   &&jPiv,
    VTAU    &&tau,
    VWORK   &&work,
    VRWORK  &&rWork)
{
    using std::min;
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type        MatrixA;
    typedef typename MatrixA::IndexType         IndexType;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();
    const IndexType k = min(m, n);

#   ifndef NDEBUG
//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(jPiv.firstIndex()==1);
    ASSERT(tau.firstIndex()==1);
    ASSERT(work.firstIndex()==1);

    ASSERT(jPiv.length()==0 || jPiv.length()==n);
    ASSERT(tau.length()==0 || tau.length()==k);
    ASSERT(work.length()>=n+1 || work.length()==IndexType(0));
    ASSERT(rWork.length()==0 || rWork.length()>=2*n);
#   endif

//
//  Resize vector output arguments and workspace if needed
//
    if (jPiv.length()==0) {
        jPiv.resize(n);
    }
    if (tau.length()==0) {
        tau.resize(k);
    }
    if (work.length()==0) {
        work.resize(qp3_wsq(A));
    }
    if (rWork.length()==0) {
        rWork.resize(2*n);
    }

//
//  Call implementation
//
    external::qp3_impl(A, jPiv, tau, work, rWork);
}

//
//-- (ge)qp3 [complex variant (with temporary workspace) ] ---------------------
//
template <typename MA, typename VJPIV, typename VTAU>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsIntegerDenseVector<VJPIV>::value
                 && IsComplexDenseVector<VTAU>::value,
         void>::Type
qp3(MA      &&A,
    VJPIV   &&jPiv,
    VTAU    &&tau)
{
    typedef typename RemoveRef<MA>::Type::Vector    Vector;

    typedef typename RemoveRef<MA>::Type::ElementType   T;
    typedef typename ComplexTrait<T>::PrimitiveType     PT;
    typedef DenseVector<Array<PT> >                     RealVector;

    Vector      work;
    RealVector  realWork;


    qp3(A, jPiv, tau, work, realWork);
}

#endif // USE_CXXLAPACK


//== workspace query ===========================================================

template <typename MA>
typename RestrictTo<IsGeMatrix<MA>::value,
         typename MA::IndexType>::Type
qp3_wsq(const MA &A)
{
    typedef typename MA::IndexType     IndexType;
//
//  Call implementation
//
    const IndexType info = LAPACK_SELECT::qp3_wsq_impl(A);

#   ifdef CHECK_CXXLAPACK
//
//  Compare generic results with results from the native implementation
//
    const IndexType _info = external::qp3_wsq_impl(A);

    ASSERT(info==_info);
#   endif
    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GE_QP3_TCC
