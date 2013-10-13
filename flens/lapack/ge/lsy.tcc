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
       SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,
      $                   WORK, LWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_GE_LSY_TCC
#define FLENS_LAPACK_GE_LSY_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- (ge)lsy [real variant] ----------------------------------------------------

template <typename MA, typename MB, typename VJPIV, typename RCOND,
          typename RANK, typename VWORK>
void
lsy_impl(GeMatrix<MA>              &A,
         GeMatrix<MB>              &B,
         DenseVector<VJPIV>        &jPiv,
         RCOND                     rCond,
         RANK                      &rank,
         DenseVector<VWORK>        &work)
{
    using std::abs;
    using flens::max;
    using flens::min;
    using LASCL::FullMatrix;
    using LASCL::UpperTriangular;

    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const ElementType Zero(0), One(1);

    const Underscore<IndexType>  _;

    const IndexType m    = A.numRows();
    const IndexType n    = A.numCols();
    const IndexType nRhs = B.numCols();
    const IndexType mn   = min(m, n);

//
//  Figure out optimal block size
//
    IndexType lWork = work.length();
    IndexType lWorkMin;
    IndexType lWorkOpt;

    if (mn==0 || nRhs==0) {
        lWorkMin = 1;
        lWorkOpt = 1;
    } else {
        IndexType nb1 = ilaenv<ElementType>(1, "GEQRF", "", m, n);
        IndexType nb2 = ilaenv<ElementType>(1, "GERQF", "", m, n);
        IndexType nb3 = ilaenv<ElementType>(1, "ORMQR", "", m, n, nRhs);
        IndexType nb4 = ilaenv<ElementType>(1, "ORMRQ", "", m, n, nRhs);
        IndexType nb = max(nb1, nb2, nb3, nb4);
        lWorkMin = mn + max(2*mn, n+1, mn+nRhs);
        lWorkOpt = max(lWorkMin, mn+2*n+nb*(n+1), 2*mn+nb*nRhs);
    }

    if (lWork==0) {
        work.resize(lWorkOpt);
        lWork = work.length();
    }
    work(1) = lWorkOpt;

    auto sMinWork = work(_(mn+1,2*mn));
    auto sMaxWork = work(_(2*mn+1,3*mn));

//
//  Quick return if possible
//
    if (mn==0 || nRhs==0) {
        rank = 0;
        return;
    }
//
//  Get machine parameters
//
    ElementType smallNum = lamch<ElementType>(SafeMin)
                         / lamch<ElementType>(Precision);
    ElementType bigNum = One / smallNum;
    labad(smallNum, bigNum);

//
//  Scale A, B if max entries outside range [SMLNUM,BIGNUM]
//
    ElementType normA = lan(MaximumNorm, A);

    IndexType scaleA = 0;

    if (normA>Zero && normA<smallNum) {
//
//      Scale matrix norm up to SMLNUM
//
        lascl(FullMatrix, 0, 0, normA, smallNum, A);
        scaleA = 1;
    } else if (normA>bigNum) {
//
//      Scale matrix norm down to BIGNUM
//
        lascl(FullMatrix, 0, 0, normA, bigNum, A);
        scaleA = 2;
    } else if (normA==Zero) {
//
//      Matrix all zero. Return zero solution.
//
        B = Zero;
        rank = 0;
        work(1) = lWorkOpt;
        return;
    }

    auto B_ = B(_(1,m),_);
    ElementType normB = lan(MaximumNorm, B_);

    IndexType scaleB = 0;

    if (normB>Zero && normB<smallNum) {
//
//      Scale matrix norm up to SMLNUM
//
        lascl(FullMatrix, 0, 0, normB, smallNum, B_);
        scaleB = 1;
    } else if (normB>bigNum) {
//
//      Scale matrix norm down to BIGNUM
//
        lascl(FullMatrix, 0, 0, normB, bigNum, B_);
        scaleB = 2;
    }
//
//  Compute QR factorization with column pivoting of A:
//      A * P = Q * R
//
    auto tau1    = work(_(1,mn));
    auto qp3Work = work(_(mn+1,lWork));
    qp3(A, jPiv, tau1, qp3Work);
//
//  workspace: MN+2*N+NB*(N+1).
//  Details of Householder rotations stored in WORK(1:MN).
//
//  Determine RANK using incremental condition estimation
//
    sMinWork(1) = One;
    sMaxWork(1) = One;

    ElementType sMax = abs(A(1,1));
    ElementType sMin = sMax;

    if (abs(A(1,1))==Zero) {
        rank = 0;
        B = Zero;
        work(1) = lWorkOpt;
        return;
    } else {
        rank = 1;
    }

    while (rank<mn) {
        IndexType i = rank+1;
        const auto sMinWork_ = sMinWork(_(1,rank));
        const auto sMaxWork_ = sMaxWork(_(1,rank));
        const auto A_i       = A(_(1,rank),i);

        ElementType sMinPr, sMaxPr, s1, s2, c1, c2;

        laic1(LAIC1::Min, sMinWork_, sMin, A_i, A(i,i), sMinPr, s1, c1);
        laic1(LAIC1::Max, sMaxWork_, sMax, A_i, A(i,i), sMaxPr, s2, c2);

        if (sMaxPr*rCond<=sMinPr) {
            for (IndexType i=1; i<=rank; ++i) {
                sMinWork(i) *= s1;
                sMaxWork(i) *= s2;
            }
            sMinWork(rank+1) = c1;
            sMaxWork(rank+1) = c2;
            sMin = sMinPr;
            sMax = sMaxPr;
            ++rank;
        } else {
            break;
        }
    }

//
//  workspace: 3*MN.
//
//  Logically partition R = [ R11 R12 ]
//                          [  0  R22 ]
//  where R11 = R(1:RANK,1:RANK)
//
//  [R11,R12] = [ T11, 0 ] * Y
//
    auto tau2 = work(_(mn+1,mn+rank));
    auto work_ = work(_(2*mn+1,lWork));
    if (rank<n) {
        auto A_ = A(_(1,rank),_);
        tzrzf(A_, tau2, work_);
    }
    auto T11 = A(_(1,rank),_(1,rank)).upper();
//
//  workspace: 2*MN.
//  Details of Householder rotations stored in WORK(MN+1:2*MN)
//
//  B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
//
    ormqr(Left, Trans, A(_,_(1,mn)), tau1, B_, work_);
//
//  workspace: 2*MN+NB*NRHS.
//
//  B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)
//
    blas::sm(Left, NoTrans, One, T11, B(_(1,rank),_));

    B(_(rank+1,n),_) = Zero;
//
//  B(1:N,1:NRHS) := Y**T * B(1:N,1:NRHS)
//
    if (rank<n) {
        ormrz(Left, Trans, n-rank, A(_(1,rank),_), tau2, B(_(1,n),_), work_);
    }
//
//  workspace: 2*MN+NRHS.
//
//  B(1:N,1:NRHS) := P * B(1:N,1:NRHS)
//
    for (IndexType j=1; j<=nRhs; ++j) {
        for (IndexType i=1; i<=n; ++i) {
            work(jPiv(i)) = B(i,j);
        }
        B(_(1,n),j) = work(_(1,n));
    }

//
//  workspace: N
//
//  Undo scaling
//
    if (scaleA==1) {
        lascl(FullMatrix, 0, 0, normA, smallNum, B(_(1,n),_));
        lascl(UpperTriangular, 0, 0, smallNum, normA, A(_(1,rank),_(1,rank)));
    } else if (scaleA==2) {
        lascl(FullMatrix, 0, 0, normA, bigNum, B(_(1,n),_));
        lascl(UpperTriangular, 0, 0, bigNum, normA, A(_(1,rank),_(1,rank)));
    }
    if (scaleB==1) {
        lascl(FullMatrix, 0, 0, smallNum, normB, B(_(1,n),_));
    } else if (scaleB==2) {
        lascl(FullMatrix, 0, 0, bigNum, normB, B(_(1,n),_));
    }

    work(1) = lWorkOpt;
}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (ge)lsy [real variant] ----------------------------------------------------

template <typename MA, typename MB, typename VJPIV, typename RCOND,
          typename RANK, typename VWORK>
void
lsy_impl(GeMatrix<MA>              &A,
         GeMatrix<MB>              &B,
         DenseVector<VJPIV>        &jPiv,
         RCOND                     rCond,
         RANK                      &rank,
         DenseVector<VWORK>        &work)
{
    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    if (work.length()==0) {
        ElementType   WORK;
        IndexType     LWORK = -1;

        cxxlapack::gelsy<IndexType>(A.numRows(),
                                    A.numCols(),
                                    B.numCols(),
                                    A.data(),
                                    A.leadingDimension(),
                                    B.data(),
                                    B.leadingDimension(),
                                    jPiv.data(),
                                    rCond,
                                    rank,
                                    &WORK,
                                    LWORK);
        work.resize(IndexType(WORK));
    }

    cxxlapack::gelsy<IndexType>(A.numRows(),
                                A.numCols(),
                                B.numCols(),
                                A.data(),
                                A.leadingDimension(),
                                B.data(),
                                B.leadingDimension(),
                                jPiv.data(),
                                rCond,
                                rank,
                                work.data(),
                                work.length());
}

//-- (ge)lsy [complex variant] -------------------------------------------------

template <typename MA, typename MB, typename VJPIV, typename RCOND,
          typename RANK, typename VWORK, typename VRWORK>
void
lsy_impl(GeMatrix<MA>              &A,
         GeMatrix<MB>              &B,
         DenseVector<VJPIV>        &jPiv,
         RCOND                     rCond,
         RANK                      &rank,
         DenseVector<VWORK>        &work,
         DenseVector<VRWORK>       &rwork)
{
    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    if (work.length()==0) {
        ElementType   WORK;
        IndexType     LWORK = -1;

        cxxlapack::gelsy<IndexType>(A.numRows(),
                                    A.numCols(),
                                    B.numCols(),
                                    A.data(),
                                    A.leadingDimension(),
                                    B.data(),
                                    B.leadingDimension(),
                                    jPiv.data(),
                                    rCond,
                                    rank,
                                    &WORK,
                                    LWORK,
                                    rwork.data());
        work.resize(IndexType(cxxblas::real(WORK)));
    }

    cxxlapack::gelsy<IndexType>(A.numRows(),
                                A.numCols(),
                                B.numCols(),
                                A.data(),
                                A.leadingDimension(),
                                B.data(),
                                B.leadingDimension(),
                                jPiv.data(),
                                rCond,
                                rank,
                                work.data(),
                                work.length(),
                                rwork.data());
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- (ge)lsy [real variant] ----------------------------------------------------

template <typename MA, typename MB, typename VJPIV, typename RCOND,
          typename VWORK>
typename RestrictTo<IsRealGeMatrix<MA>::value
                 && IsRealGeMatrix<MB>::value
                 && IsIntegerDenseVector<VJPIV>::value
                 && IsReal<RCOND>::value
                 && IsRealDenseVector<VWORK>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
lsy(MA           &&A,
    MB           &&B,
    VJPIV        &&jPiv,
    RCOND        rCond,
    VWORK        &&work)
{
    using flens::max;
    using std::min;

    LAPACK_DEBUG_OUT("(ge)lsy [real]");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MB>::Type    MatrixB;
    typedef typename RemoveRef<VJPIV>::Type VectorJPiv;
    typedef typename RemoveRef<VWORK>::Type VectorWork;
#   endif

//
//  Test the input parameters
//
    const IndexType n = A.numCols();

#   ifndef NDEBUG
    const IndexType m = A.numRows();
    const IndexType nRhs = B.numCols();

    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(jPiv.firstIndex()==1);
    ASSERT(work.firstIndex()==1);
    ASSERT(B.numRows()==max(m,n));
    ASSERT(jPiv.length()==0 || jPiv.length()==n);

    if (work.length()>0) {
        const IndexType mn = min(m, n);
        const IndexType lWorkMin = mn + max(2*mn, n + 1, mn + nRhs);
        ASSERT(work.length()>=lWorkMin);
    }
#   endif

    if (jPiv.length()==0) {
        jPiv.resize(n, jPiv.firstIndex(), IndexType(0));
    }

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView        A_org      = A;
    typename MatrixB::NoView        B_org      = B;
    typename VectorJPiv::NoView     jPiv_org   = jPiv;
    typename VectorWork::NoView     work_org   = work;
#   endif

//
//  Call implementation
//
    IndexType  rank;
    LAPACK_SELECT::lsy_impl(A, B, jPiv, rCond, rank, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename MatrixA::NoView        A_generic    = A;
    typename MatrixB::NoView        B_generic    = B;
    typename VectorJPiv::NoView     jPiv_generic = jPiv;
    IndexType                       rank_generic = rank;
    typename VectorWork::NoView     work_generic = work;

//
//  restore output arguments
//
    A    = A_org;
    B    = B_org;
    jPiv = jPiv_org;
    rank = 0;
    work = work_org;

//
//  Compare generic results with results from the native implementation
//
    external::lsy_impl(A, B, jPiv, rCond, rank, work);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }
    if (! isIdentical(B_generic, B, "B_generic", "B")) {
        std::cerr << "CXXLAPACK: B_generic = " << B_generic << std::endl;
        std::cerr << "F77LAPACK: B = " << B << std::endl;
        failed = true;
    }
    if (! isIdentical(jPiv_generic, jPiv, "jPiv_generic", "jPiv")) {
        std::cerr << "CXXLAPACK: jPiv_generic = " << jPiv_generic << std::endl;
        std::cerr << "F77LAPACK: jPiv = " << jPiv << std::endl;
        failed = true;
    }
    if (! isIdentical(rank_generic, rank, "rank_generic", "rank")) {
        std::cerr << "CXXLAPACK: rank_generic = " << rank_generic << std::endl;
        std::cerr << "F77LAPACK: rank = " << rank << std::endl;
        failed = true;
    }
    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "A.numRows() = " << A.numRows() << std::endl;
        std::cerr << "A.numCols() = " << A.numCols() << std::endl;
        std::cerr << "A = " << A << std::endl;
        std::cerr << "rank_generic = " << rank_generic << std::endl;
        std::cerr << "rank = " << rank << std::endl;
        std::cerr << "error in: lsy.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: lsy.tcc" << std::endl;
    }
#   endif

    return rank;
}

//-- (ge)lsy [complex variant] -------------------------------------------------

#ifdef USE_CXXLAPACK

template <typename MA, typename MB, typename VJPIV, typename RCOND,
          typename VWORK, typename VRWORK>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexGeMatrix<MB>::value
                 && IsIntegerDenseVector<VJPIV>::value
                 && IsReal<RCOND>::value
                 && IsComplexDenseVector<VWORK>::value
                 && IsRealDenseVector<VRWORK>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
lsy(MA           &&A,
    MB           &&B,
    VJPIV        &&jPiv,
    RCOND        rCond,
    VWORK        &&work,
    VRWORK       &&rwork)
{
    using flens::max;
    using flens::min;

    LAPACK_DEBUG_OUT("(ge)lsy [complex]");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

//
//  Test the input parameters
//
    const IndexType n = A.numCols();

#   ifndef NDEBUG
    const IndexType m = A.numRows();
    const IndexType nRhs = B.numCols();

    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(jPiv.firstIndex()==1);
    ASSERT(work.firstIndex()==1);
    ASSERT(B.numRows()==max(m,n));
    ASSERT(jPiv.length()==0 || jPiv.length()==n);

    if (work.length()>0) {
        const IndexType mn = min(m, n);
        const IndexType lWorkMin = mn + max(2*mn, n + 1, mn + nRhs);
        ASSERT(work.length()>=lWorkMin);
    }

    ASSERT(rwork.length()==0 || rwork.length()==2*n);
#   endif

    if (jPiv.length()==0) {
        jPiv.resize(n, jPiv.firstIndex(), IndexType(0));
    }

    if (rwork.length()==0) {
        rwork.resize(2*n);
    }

//
//  Call implementation
//
    IndexType  rank;
    external::lsy_impl(A, B, jPiv, rCond, rank, work, rwork);

    return rank;
}

#endif // USE_CXXLAPACK


//-- (ge)lsy [real variant with temporary workspace] ---------------------------

template <typename MA, typename MB, typename VJPIV, typename RCOND>
typename RestrictTo<IsRealGeMatrix<MA>::value
                 && IsRealGeMatrix<MB>::value
                 && IsIntegerDenseVector<VJPIV>::value
                 && IsReal<RCOND>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
lsy(MA           &&A,
    MB           &&B,
    VJPIV        &&jPiv,
    RCOND        rCond)
{
    typedef typename RemoveRef<MA>::Type::Vector WorkVector;

    WorkVector  work;
    return lsy(A, B, jPiv, rCond);
}


//-- (ge)lsy [complex variant with temporary workspace] ------------------------

#ifdef USE_CXXLAPACK

template <typename MA, typename MB, typename VJPIV, typename RCOND>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexGeMatrix<MB>::value
                 && IsIntegerDenseVector<VJPIV>::value
                 && IsReal<RCOND>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
lsy(MA           &&A,
    MB           &&B,
    VJPIV        &&jPiv,
    RCOND        rCond)
{
    typedef typename RemoveRef<MA>::Type::Vector        WorkVector;
    typedef typename RemoveRef<MA>::Type::ElementType   T;
    typedef typename ComplexTrait<T>::PrimitiveType     PT;
    typedef DenseVector<Array<PT> >                     RealWorkVector;

    WorkVector      work;
    RealWorkVector  rwork;
    return lsy(A, B, jPiv, rCond, work, rwork);
}

#endif // USE_CXXLAPACK


//== (ge)lsy variant if B is vector ============================================


//-- (ge)lsy [real variant] ----------------------------------------------------

template <typename MA, typename VB, typename VJPIV, typename RCOND,
          typename VWORK>
typename RestrictTo<IsRealGeMatrix<MA>::value
                 && IsRealDenseVector<VB>::value
                 && IsIntegerDenseVector<VJPIV>::value
                 && IsReal<RCOND>::value
                 && IsRealDenseVector<VWORK>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
lsy(MA           &&A,
    VB           &&b,
    VJPIV        &&jPiv,
    RCOND        rCond,
    VWORK        &&work)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VB>::Type    VectorB;

    typedef typename VectorB::ElementType  ElementType;
    typedef typename VectorB::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = MatrixA::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b, n);
    return lsy(A, B, jPiv, rCond, work);
}


//-- (ge)lsy [complex variant] -------------------------------------------------

#ifdef USE_CXXLAPACK

template <typename MA, typename VB, typename VJPIV, typename RCOND,
          typename VWORK, typename VRWORK>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexDenseVector<VB>::value
                 && IsIntegerDenseVector<VJPIV>::value
                 && IsReal<RCOND>::value
                 && IsComplexDenseVector<VWORK>::value
                 && IsRealDenseVector<VRWORK>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
lsy(MA           &&A,
    VB           &&b,
    VJPIV        &&jPiv,
    RCOND        rCond,
    VWORK        &&work,
    VRWORK       &&rwork)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VB>::Type    VectorB;

    typedef typename VectorB::ElementType  ElementType;
    typedef typename VectorB::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = MatrixA::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b, n);
    return lsy(A, B, jPiv, rCond, work, rwork);
}

#endif // USE_CXXLAPACK

//-- (ge)lsy [real variant with temporary workspace] ---------------------------

template <typename MA, typename VB, typename VJPIV, typename RCOND>
typename RestrictTo<IsRealGeMatrix<MA>::value
                 && IsRealDenseVector<VB>::value
                 && IsIntegerDenseVector<VJPIV>::value
                 && IsReal<RCOND>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
lsy(MA           &&A,
    VB           &&b,
    VJPIV        &&jPiv,
    RCOND        rCond)
{
    typedef typename RemoveRef<MA>::Type::Vector WorkVector;

    WorkVector  work;
    return lsy(A, b, jPiv, rCond, work);
}


//-- (ge)lsy [complex variant with temporary workspace] ------------------------

#ifdef USE_CXXLAPACK

template <typename MA, typename VB, typename VJPIV, typename RCOND>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexDenseVector<VB>::value
                 && IsIntegerDenseVector<VJPIV>::value
                 && IsReal<RCOND>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
lsy(MA           &&A,
    VB           &&b,
    VJPIV        &&jPiv,
    RCOND        rCond)
{
    typedef typename RemoveRef<MA>::Type::Vector        WorkVector;
    typedef typename RemoveRef<MA>::Type::ElementType   T;
    typedef typename ComplexTrait<T>::PrimitiveType     PT;
    typedef DenseVector<Array<PT> >                     RealWorkVector;

    WorkVector      work;
    RealWorkVector  rwork;
    return lsy(A, b, jPiv, rCond, work, rwork);
}

#endif // USE_CXXLAPACK



} } // namespace lapack, flens

#endif // FLENS_LAPACK_GE_LSY_TCC
