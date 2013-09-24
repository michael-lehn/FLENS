/*
 *   Copyright (c) 2012, Michael Lehn, Klaus Pototzky
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
      SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
    $                    WORK, LWORK, INFO )
      SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
    $                    WORK, LWORK, RWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_GE_SVD_TCC
#define FLENS_LAPACK_GE_SVD_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (ge)svd [real variant] ----------------------------------------

template <typename MA>
typename RestrictTo<IsRealGeMatrix<MA>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
svd_wsq_impl(SVD::Job           jobU,
             SVD::Job           jobVT,
             GeMatrix<MA>       &A)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;
    typedef typename GeMatrix<MA>::ElementType  T;

    T               WORK, DUMMY;
    const IndexType LWORK   = -1;
    const IndexType ONE     =  1;

    cxxlapack::gesvd<IndexType>(getF77Char(jobU), getF77Char(jobVT),
                                A.numRows(),
                                A.numCols(),
                                A.data(),
                                A.leadingDimension(),
                                &DUMMY,
                                &DUMMY,
                                ONE,
                                &DUMMY,
                                ONE,
                                &WORK,
                                LWORK);

    return WORK;
}

template <typename MA, typename VS, typename MU, typename MVT, typename VWORK>
void
svd_impl(SVD::Job               jobU,
         SVD::Job               jobVT,
         GeMatrix<MA>           &A,
         DenseVector<VS>        &s,
         GeMatrix<MU>           &U,
         GeMatrix<MVT>          &VT,
         DenseVector<VWORK>     &work)
{
    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    if (work.length()==0) {
        ElementType     WORK;
        IndexType       LWORK = -1;

        cxxlapack::gesvd<IndexType>(getF77Char(jobU), getF77Char(jobVT),
                                    A.numRows(),
                                    A.numCols(),
                                    A.data(),
                                    A.leadingDimension(),
                                    s.data(),
                                    U.data(),
                                    U.leadingDimension(),
                                    VT.data(),
                                    VT.leadingDimension(),
                                    &WORK,
                                    LWORK);
        work.resize(cxxblas::real(WORK));
    }
    cxxlapack::gesvd<IndexType>(getF77Char(jobU), getF77Char(jobVT),
                                A.numRows(),
                                A.numCols(),
                                A.data(),
                                A.leadingDimension(),
                                s.data(),
                                U.data(),
                                U.leadingDimension(),
                                VT.data(),
                                VT.leadingDimension(),
                                work.data(),
                                work.length());
}

//-- (ge)svd [complex variant] ----------------------------------------

template <typename MA>
typename RestrictTo<IsComplexGeMatrix<MA>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
svd_wsq_impl(SVD::Job           jobU,
             SVD::Job           jobVT,
             GeMatrix<MA>       &A)
{
    typedef typename GeMatrix<MA>::IndexType        IndexType;
    typedef typename GeMatrix<MA>::ElementType      T;
    typedef typename ComplexTrait<T>::PrimitiveType PT;

    T               WORK, DUMMY;
    PT              RWORK, RDUMMY;
    const IndexType LWORK   = -1;
    const IndexType ONE     =  1;

    cxxlapack::gesvd<IndexType>(getF77Char(jobU), getF77Char(jobVT),
                                A.numRows(),
                                A.numCols(),
                                A.data(),
                                A.leadingDimension(),
                                &RDUMMY,
                                &DUMMY,
                                ONE,
                                &DUMMY,
                                ONE,
                                &WORK,
                                LWORK,
                                &RWORK);

    return cxxblas::real(WORK);
}

template <typename MA, typename VS, typename MU, typename MVT,
          typename VWORK, typename VRWORK>
void
svd_impl(SVD::Job               jobU,
         SVD::Job               jobVT,
         GeMatrix<MA>           &A,
         DenseVector<VS>        &s,
         GeMatrix<MU>           &U,
         GeMatrix<MVT>          &VT,
         DenseVector<VWORK>     &work,
         DenseVector<VRWORK>    &rwork)
{
    using std::min;

    typedef typename GeMatrix<MA>::ElementType      ElementType;
    typedef typename GeMatrix<MA>::IndexType        IndexType;

    if (work.length()==0) {
        ElementType     WORK;
        IndexType       LWORK = -1;

        cxxlapack::gesvd<IndexType>(getF77Char(jobU), getF77Char(jobVT),
                                    A.numRows(),
                                    A.numCols(),
                                    A.data(),
                                    A.leadingDimension(),
                                    s.data(),
                                    U.data(),
                                    U.leadingDimension(),
                                    VT.data(),
                                    VT.leadingDimension(),
                                    &WORK,
                                    LWORK,
                                    rwork.data());
        work.resize(cxxblas::real(WORK));
    }

    if (rwork.length()==0) {
        rwork.resize(5*min(A.numRows(), A.numCols()));
    }

    typedef typename GeMatrix<MA>::IndexType  IndexType;

    cxxlapack::gesvd<IndexType>(getF77Char(jobU), getF77Char(jobVT),
                                A.numRows(),
                                A.numCols(),
                                A.data(),
                                A.leadingDimension(),
                                s.data(),
                                U.data(),
                                U.leadingDimension(),
                                VT.data(),
                                VT.leadingDimension(),
                                work.data(),
                                work.length(),
                                rwork.data());
}

} // namespace external

//== public interface ==========================================================

//-- (ge)svd [real variant] ----------------------------------------------------

template <typename MA, typename VS, typename MU, typename MVT, typename VWORK>
typename RestrictTo<IsRealGeMatrix<MA>::value
                 && IsRealDenseVector<VS>::value
                 && IsRealGeMatrix<MU>::value
                 && IsRealGeMatrix<MVT>::value
                 && IsRealDenseVector<VWORK>::value,
          void>::Type
svd(SVD::Job    jobU,
    SVD::Job    jobVT,
    MA          &&A,
    VS          &&s,
    MU          &&U,
    MVT         &&VT,
    VWORK       &&work)
{
    using std::min;

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type      MatrixA;
    typedef typename MatrixA::IndexType       IndexType;


#   ifndef NDEBUG
//
//  Test the input parameters
//
    ASSERT((jobU != SVD::Overwrite) || (jobVT != SVD::Overwrite));
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(s.firstIndex()==1);
    ASSERT(work.firstIndex()==1);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    ASSERT(s.length()==std::min(m,n));
#   endif

//
//  Call implementation
//
    external::svd_impl(jobU, jobVT, A, s, U, VT, work);
}


//-- (ge)svd [complex variant] -------------------------------------------------


template <typename MA, typename VS, typename MU, typename MVT, typename VWORK,
          typename VRWORK>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsRealDenseVector<VS>::value
                 && IsComplexGeMatrix<MU>::value
                 && IsComplexGeMatrix<MVT>::value
                 && IsComplexDenseVector<VWORK>::value
                 && IsRealDenseVector<VRWORK>::value,
          void>::Type
svd(SVD::Job    jobU,
    SVD::Job    jobVT,
    MA          &&A,
    VS          &&s,
    MU          &&U,
    MVT         &&VT,
    VWORK       &&work,
    VRWORK      &&rwork)
{
    using std::min;

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type      MatrixA;
    typedef typename MatrixA::IndexType       IndexType;

#   ifndef NDEBUG
//
//  Test the input parameters
//
    ASSERT((jobU != SVD::Overwrite) || (jobVT != SVD::Overwrite));
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(s.firstIndex()==1);
    ASSERT(work.firstIndex()==1);
    ASSERT(rwork.firstIndex()==1);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    ASSERT(s.length()==std::min(m,n));
#   endif

//
//  Call implementation
//
    external::svd_impl(jobU, jobVT, A, s, U, VT, work, rwork);
}


//
//  Real variant with temporary workspace
//
template <typename MA, typename VS, typename MU, typename MVT>
typename RestrictTo<IsRealGeMatrix<MA>::value
                 && IsRealDenseVector<VS>::value
                 && IsRealGeMatrix<MU>::value
                 && IsRealGeMatrix<MVT>::value,
          void>::Type
svd(SVD::Job    jobU,
    SVD::Job    jobVT,
    MA          &&A,
    VS          &&s,
    MU          &&U,
    MVT         &&VT)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type::Vector WorkVector;
    WorkVector work;
    svd(jobU, jobVT, A, s, U, VT, work);
}


//
//  Complex variant with temporary workspace
//
template <typename MA, typename VS, typename MU, typename MVT>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsRealDenseVector<VS>::value
                 && IsComplexGeMatrix<MU>::value
                 && IsComplexGeMatrix<MVT>::value,
          void>::Type
svd(SVD::Job    jobU,
    SVD::Job    jobVT,
    MA          &&A,
    VS          &&s,
    MU          &&U,
    MVT         &&VT)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type::Vector        WorkVector;
    typedef typename RemoveRef<MA>::Type::ElementType   T;
    typedef typename ComplexTrait<T>::PrimitiveType     PT;
    typedef DenseVector<Array<PT> >                     RealWorkVector;

    WorkVector      work;
    RealWorkVector  rwork;

    svd(jobU, jobVT, A, s, U, VT, work, rwork);
}


//-- svd_wsq [worksize query] ------------------------------------------------
template <typename MA, typename VS, typename MU, typename MVT, typename VWORK>
typename RestrictTo<IsGeMatrix<MA>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
svd_wsq(SVD::Job    jobU,
        SVD::Job    jobVT,
        MA          &&A)
{
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT((jobU != SVD::Overwrite) || (jobVT != SVD::Overwrite));
#   endif

//
//  Call implementation
//
    const IndexType info = external::svd_wsq_impl(jobU, jobVT, A);

    return info;
}

#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GE_SVD_TCC
