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
       SUBROUTINE DGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,
      $                   X, LDX, FERR, BERR, WORK, IWORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_GESV_RFS_TCC
#define FLENS_LAPACK_GESV_RFS_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
template <typename MA, typename MAF, typename VPIV, typename MB, typename MX,
          typename VFERR, typename VBERR, typename VWORK, typename VIWORK>
void
rfs_generic(Transpose               trans,
            const GeMatrix<MA>      &A,
            const GeMatrix<MAF>     &AF,
            const DenseVector<VPIV> &piv,
            const GeMatrix<MB>      &B,
            GeMatrix<MX>            &X,
            DenseVector<VFERR>      &fErr,
            DenseVector<VBERR>      &bErr,
            DenseVector<VWORK>      &work,
            DenseVector<VIWORK>     &iwork)
{
    using std::abs;

    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const IndexType   itMax = 5;
    const ElementType Zero(0), One(1), Two(2), Three(3);

    const Underscore<IndexType> _;

    const IndexType n    = B.numRows();
    const IndexType nRhs = B.numCols();
//
//  Local Arrays
//
    IndexType iSaveData[3] = {0, 0, 0};
    DenseVectorView<IndexType>
        iSave = typename DenseVectorView<IndexType>::Engine(3, iSaveData, 1);

//
//  Quick return if possible
//
    if (n==0 || nRhs==0) {
        fErr = Zero;
        bErr = Zero;
        return;
    }

    const Transpose transT = (trans==NoTrans) ? Trans : NoTrans;
//
//  NZ = maximum number of nonzero elements in each row of A, plus 1
//
    const IndexType nz = n + 1;
    const ElementType eps = lamch<ElementType>(Eps);
    const ElementType safeMin = lamch<ElementType>(SafeMin);
    const ElementType safe1 = nz * safeMin;
    const ElementType safe2 = safe1 / eps;

    auto work1 = work(_(1,n));
    auto work2 = work(_(n+1,2*n));
    auto work3 = work(_(2*n+1,3*n));
//
//  Do for each right hand side
//

    for (IndexType j=1; j<=nRhs; ++j) {
        IndexType   count = 1;
        ElementType lastRes = Three;

        RETRY:
//
//      Loop until stopping criterion is satisfied.
//
//      Compute residual R = B - op(A) * X,
//      where op(A) = A, A**T, or A**H, depending on TRANS.
//
        work2 = B(_,j);
        blas::mv(trans, -One, A, X(_,j), One, work2);
//
//      Compute componentwise relative backward error from formula
//
//      max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )
//
//      where abs(Z) is the componentwise absolute value of the matrix
//      or vector Z.  If the i-th component of the denominator is less
//      than SAFE2, then SAFE1 is added to the i-th components of the
//      numerator and denominator before dividing.
//
        for (IndexType i=1; i<=n; ++i) {
            work(i) = abs(B(i,j));
        }
//
//      Compute abs(op(A))*abs(X) + abs(B).
//
        if (trans==NoTrans) {
            for (IndexType k=1; k<=n; ++k) {
                const ElementType xk = abs(X(k,j));
                for (IndexType i=1; i<=n; ++i) {
                    work(i) += abs(A(i,k)) * xk;
                }
            }
        } else {
            for (IndexType k=1; k<=n; ++k) {
                ElementType s = Zero;
                for (IndexType i=1; i<=n; ++i) {
                    s += abs(A(i,k)) * abs(X(i,j));
                }
                work1(k) += s;
            }
        }

        ElementType s = Zero;
        for (IndexType i=1; i<=n; ++i) {
            if (work1(i)>safe2) {
                s = max(s, abs(work2(i))/work1(i));
            } else {
                s = max(s, (abs(work2(i))+safe1)/(work1(i)+safe1));
            }
        }
        bErr(j) = s;
//
//      Test stopping criterion. Continue iterating if
//         1) The residual BERR(J) is larger than machine epsilon, and
//         2) BERR(J) decreased by at least a factor of 2 during the
//            last iteration, and
//         3) At most ITMAX iterations tried.
//

        if (bErr(j)>eps && Two*bErr(j)<=lastRes && count<=itMax) {
//
//          Update solution and try again.
//
            trs(trans, AF, piv, work2);
            X(_,j) += work2;
            lastRes = bErr(j);
            ++count;
            goto RETRY;
        }
//
//      Bound error from formula
//
//      norm(X - XTRUE) / norm(X) .le. FERR =
//      norm( abs(inv(op(A)))*
//         ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)
//
//      where
//        norm(Z) is the magnitude of the largest component of Z
//        inv(op(A)) is the inverse of op(A)
//        abs(Z) is the componentwise absolute value of the matrix or
//           vector Z
//        NZ is the maximum number of nonzeros in any row of A, plus 1
//        EPS is machine epsilon
//
//      The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
//      is incremented by SAFE1 if the i-th component of
//      abs(op(A))*abs(X) + abs(B) is less than SAFE2.
//
//      Use DLACN2 to estimate the infinity-norm of the matrix
//         inv(op(A)) * diag(W),
//      where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))
//
        for (IndexType i=1; i<=n; ++i) {
            if (work(i)>safe2) {
                work(i) = abs(work2(i)) + nz*eps*work1(i);
            } else {
                work(i) = abs(work2(i)) + nz*eps*work1(i) + safe1;
            }
        }

        IndexType kase = 0;
        while (true) {
            lacn2(work3, work2, iwork, fErr(j), kase, iSave);
            if (kase==0) {
                break;
            }
            if (kase==1) {
//
//              Multiply by diag(W)*inv(op(A)**T).
//
                trs(transT, AF, piv, work2);
                for (IndexType i=1; i<=n; ++i) {
                    work2(i) *= work1(i);
                }
            } else {
//
//              Multiply by inv(op(A))*diag(W).
//
                for (IndexType i=1; i<=n; ++i) {
                    work2(i) *= work1(i);
                }
                trs(trans, AF, piv, work2);
            }
        }
//
//      Normalize error.
//
        lastRes = Zero;
        for (IndexType i=1; i<=n; ++i) {
            lastRes = max(lastRes, abs(X(i,j)));
        }
        if (lastRes!=Zero) {
            fErr(j) /= lastRes;
        }

    }
}

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename MAF, typename VPIV, typename MB, typename MX,
          typename VFERR, typename VBERR, typename VWORK, typename VIWORK>
void
rfs(Transpose               trans,
    const GeMatrix<MA>      &A,
    const GeMatrix<MAF>     &AF,
    const DenseVector<VPIV> &piv,
    const GeMatrix<MB>      &B,
    GeMatrix<MX>            &X,
    DenseVector<VFERR>      &fErr,
    DenseVector<VBERR>      &bErr,
    DenseVector<VWORK>      &work,
    DenseVector<VIWORK>     &iwork)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::gerfs<IndexType>(getF77Char(trans),
                                                 B.numRows(),
                                                 B.numCols(),
                                                 A.data(),
                                                 A.leadingDimension(),
                                                 AF.data(),
                                                 AF.leadingDimension(),
                                                 piv.data(),
                                                 B.data(),
                                                 B.leadingDimension(),
                                                 X.data(),
                                                 X.leadingDimension(),
                                                 fErr.data(),
                                                 bErr.data(),
                                                 work.data(),
                                                 iwork.data());
    ASSERT(info==0);
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
template <typename MA, typename MAF, typename VPIV, typename MB, typename MX,
          typename VFERR, typename VBERR, typename VWORK, typename VIWORK>
void
rfs(Transpose               trans,
    const GeMatrix<MA>      &A,
    const GeMatrix<MAF>     &AF,
    const DenseVector<VPIV> &piv,
    const GeMatrix<MB>      &B,
    GeMatrix<MX>            &X,
    DenseVector<VFERR>      &fErr,
    DenseVector<VBERR>      &bErr,
    DenseVector<VWORK>      &work,
    DenseVector<VIWORK>     &iwork)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());

    const IndexType n = A.numRows();

    ASSERT(AF.firstRow()==1);
    ASSERT(AF.firstCol()==1);
    ASSERT(AF.numRows()==n);
    ASSERT(AF.numCols()==n);

    ASSERT(piv.firstIndex()==1);
    ASSERT(piv.length()==n);

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==n);

    const IndexType nRhs = B.numCols();

    ASSERT(X.firstRow()==1);
    ASSERT(X.firstCol()==1);
    ASSERT(X.numRows()==n);
    ASSERT(X.numCols()==nRhs);

    ASSERT(fErr.firstIndex()==1);
    ASSERT(fErr.length()==nRhs);

    ASSERT(bErr.firstIndex()==1);
    ASSERT(bErr.length()==nRhs);

    ASSERT(work.firstIndex()==1);
    ASSERT(work.length()==3*n);

    ASSERT(iwork.firstIndex()==1);
    ASSERT(iwork.length()==n);
#   endif

//
//  Make copies of output arguments
//
    typename GeMatrix<MX>::NoView        X_org     = X;
    typename DenseVector<VFERR>::NoView  fErr_org  = fErr;
    typename DenseVector<VBERR>::NoView  bErr_org  = bErr;
    typename DenseVector<VWORK>::NoView  work_org  = work;
    typename DenseVector<VIWORK>::NoView iwork_org = iwork;
//
//  Call implementation
//
    rfs_generic(trans, A, AF, piv, B, X, fErr, bErr, work, iwork);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename GeMatrix<MX>::NoView        X_generic     = X;
    typename DenseVector<VFERR>::NoView  fErr_generic  = fErr;
    typename DenseVector<VBERR>::NoView  bErr_generic  = bErr;
    typename DenseVector<VWORK>::NoView  work_generic  = work;
    typename DenseVector<VIWORK>::NoView iwork_generic = iwork;

    X     = X_org;
    fErr  = fErr_org;
    bErr  = bErr_org;
    work  = work_org;
    iwork = iwork_org;

    external::rfs(trans, A, AF, piv, B, X, fErr, bErr, work, iwork);

    bool failed = false;
    if (! isIdentical(X_generic, X, "X_generic", "X")) {
        std::cerr << "CXXLAPACK: X_generic = " << X_generic << std::endl;
        std::cerr << "F77LAPACK: X = " << X << std::endl;
        failed = true;
    }

    if (! isIdentical(fErr_generic, fErr, "fErr_generic", "fErr")) {
        std::cerr << "CXXLAPACK: fErr_generic = " << fErr_generic << std::endl;
        std::cerr << "F77LAPACK: fErr = " << fErr << std::endl;
        failed = true;
    }

    if (! isIdentical(bErr_generic, bErr, "bErr_generic", "bErr")) {
        std::cerr << "CXXLAPACK: bErr_generic = " << bErr_generic << std::endl;
        std::cerr << "F77LAPACK: bErr = " << bErr << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (! isIdentical(iwork_generic, iwork, "iwork_generic", "iwork")) {
        std::cerr << "CXXLAPACK: iwork_generic = "
                  << iwork_generic << std::endl;
        std::cerr << "F77LAPACK: iwork = " << iwork << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename MAF, typename VPIV, typename MB, typename MX,
          typename VFERR, typename VBERR, typename VWORK, typename VIWORK>
void
rfs(Transpose    trans,
    const MA     &A,
    const MAF    &AF,
    const VPIV   &piv,
    const MB     &B,
    MX           &&X,
    VFERR        &&fErr,
    VBERR        &&bErr,
    VWORK        &&work,
    VIWORK       &&iwork)
{
    CHECKPOINT_ENTER;
    rfs(trans, A, AF, piv, B, X, fErr, bErr, work, iwork);
    CHECKPOINT_LEAVE;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GESV_RFS_TCC
