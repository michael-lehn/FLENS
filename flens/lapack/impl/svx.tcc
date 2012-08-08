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
       SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
      $                   EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
      $                   WORK, IWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_SVX_TCC
#define FLENS_LAPACK_IMPL_SVX_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
namespace generic {

template <typename MA, typename MAF, typename VPIV, typename VR, typename VC,
          typename MB, typename MX, typename RCOND, typename FERR,
          typename BERR, typename VWORK, typename VIWORK>
typename GeMatrix<MA>::IndexType
svx_impl(SVX::Fact           fact,
         Transpose           trans,
         GeMatrix<MA>        &A,
         GeMatrix<MAF>       &AF,
         DenseVector<VPIV>   &piv,
         SVX::Equilibration  equed,
         DenseVector<VR>     &r,
         DenseVector<VC>     &c,
         GeMatrix<MB>        &B,
         GeMatrix<MX>        &X,
         RCOND               &rCond,
         DenseVector<FERR>   &fErr,
         DenseVector<BERR>   &bErr,
         DenseVector<VWORK>  &work,
         DenseVector<VIWORK> &iwork)
{
    using std::max;
    using std::min;
    using namespace SVX;

    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const ElementType Zero(0), One(1);

    const Underscore<IndexType> _;

    const IndexType n    = A.numRows();
    const IndexType nRhs = B.numCols();

    IndexType info = 0;

    const ElementType smallNum = lamch<ElementType>(SafeMin);
    const ElementType bigNum = One / smallNum;

    bool rowEqu, colEqu;
    ElementType rPivGrowth;

    if (fact==NotFactored || fact==Equilibrate) {
        equed = None;
        rowEqu = false;
        colEqu = false;
    } else {
        rowEqu = (equed==Row || equed==Both);
        colEqu = (equed==Column || equed==Both);
    }
//
//  Test the input parameters.
//
    ElementType rowCond, colCond;

    if (rowEqu) {
        ElementType  rcMin = bigNum;
        ElementType  rcMax = Zero;
        for (IndexType j=1; j<=n; ++j) {
            rcMin = min(rcMin, r(j));
            rcMax = max(rcMax, r(j));
        }
        ASSERT(rcMin>Zero);
        if (n>0) {
            rowCond = max(rcMin,smallNum) / min(rcMax,bigNum);
        } else {
            rowCond = One;
        }
    }
    if (colEqu) {
        ElementType  rcMin = bigNum;
        ElementType  rcMax = Zero;
        for (IndexType j=1; j<=n; ++j) {
            rcMin = min(rcMin, c(j));
            rcMax = max(rcMax, c(j));
        }
        ASSERT(rcMin>Zero);
        if (n>0) {
            colCond = max(rcMin,smallNum) / min(rcMax,bigNum);
        } else {
            colCond = One;
        }
    }

    if (fact==Equilibrate) {
        ElementType amax;
//
//      Compute row and column scalings to equilibrate the matrix A.
//
        if (equ(A, r, c, rowCond, colCond, amax)==0) {
//
//          Equilibrate the matrix.
//
            equed = Equilibration(laq(A, r, c, rowCond, colCond, amax));
            rowEqu = (equed==Row || equed==Both);
            colEqu = (equed==Column || equed==Both);
        }
    }
//
//  Scale the right hand side.
//
    if (trans==NoTrans) {
        if (rowEqu) {
            for (IndexType j=1; j<=nRhs; ++j) {
                for (IndexType i=1; i<=n; ++i) {
                    B(i,j) *= r(i);
                }
            }
        }
    } else if (colEqu) {
        for (IndexType j=1; j<=nRhs; ++j) {
            for (IndexType i=1; i<=n; ++i) {
                B(i,j) *= c(i);
            }
        }
    }

    if ((fact==NotFactored) || (fact==Equilibrate)) {
//
//      Compute the LU factorization of A.
//
        AF = A;
        info = trf(AF, piv);
//
//      Return if INFO is non-zero.
//
        if (info>0) {
//
//          Compute the reciprocal pivot growth factor of the
//          leading rank-deficient INFO columns of A.
//
            const auto range = _(1,info);
            rPivGrowth = lan(MaximumNorm, AF(range,range).upper());
            if (rPivGrowth==Zero) {
                rPivGrowth = One;
            } else {
                rPivGrowth = lan(MaximumNorm, A) / rPivGrowth;
            }
            work(1) = rPivGrowth;
            rCond = Zero;
            return info;
        }
    }
//
//  Compute the norm of the matrix A and the
//  reciprocal pivot growth factor RPVGRW.
//
    const Norm norm = (trans==NoTrans) ? OneNorm : InfinityNorm;
    const ElementType normA = lan(norm, A, work);
    rPivGrowth = lan(MaximumNorm, AF.upper());
    if (rPivGrowth==Zero) {
        rPivGrowth = One;
    } else {
        rPivGrowth = lan(MaximumNorm, A) / rPivGrowth;
    }
//
//  Compute the reciprocal of the condition number of A.
//
    con(norm, AF, normA, rCond, work, iwork);
//
//  Compute the solution matrix X.
//
    X = B;
    trs(trans, AF, piv, X);
//
//  Use iterative refinement to improve the computed solution and
//  compute error bounds and backward error estimates for it.
//
    rfs(trans, A, AF, piv, B, X, fErr, bErr, work(_(1,3*n)), iwork);
//
//  Transform the solution matrix X to a solution of the original
//  system.
//
    if (trans==NoTrans) {
        if (colEqu) {
            for (IndexType j=1; j<=nRhs; ++j) {
                for (IndexType i=1; i<=n; ++i) {
                    X(i,j) *= c(i);
                }
            }
            for (IndexType j=1; j<=nRhs; ++j) {
                fErr(j) /= colCond;
            }
        }
    } else if (rowEqu) {
        for (IndexType j=1; j<=nRhs; ++j) {
            for (IndexType i=1; i<=n; ++i) {
                X(i,j) *= r(i);
            }
        }
        for (IndexType j=1; j<=nRhs; ++j) {
            fErr(j) /= rowCond;
        }
    }

    work(1) = rPivGrowth;
//
//  Set INFO = N+1 if the matrix is singular to working precision.
//
    const ElementType eps = lamch<ElementType>(Eps);
    if (rCond<eps) {
        info = n+1;
    }
    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename MAF, typename VPIV, typename VR, typename VC,
          typename MB, typename MX, typename RCOND, typename FERR,
          typename BERR, typename VWORK, typename VIWORK>
typename GeMatrix<MA>::IndexType
svx_impl(SVX::Fact           fact,
         Transpose           trans,
         GeMatrix<MA>        &A,
         GeMatrix<MAF>       &AF,
         DenseVector<VPIV>   &piv,
         SVX::Equilibration  equed,
         DenseVector<VR>     &r,
         DenseVector<VC>     &c,
         GeMatrix<MB>        &B,
         GeMatrix<MX>        &X,
         RCOND               &rCond,
         DenseVector<FERR>   &fErr,
         DenseVector<BERR>   &bErr,
         DenseVector<VWORK>  &work,
         DenseVector<VIWORK> &iwork)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info;
    info = cxxlapack::gesvx<IndexType>(getF77Char(fact),
                                       getF77Char(trans),
                                       A.numRows(),
                                       B.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       AF.data(),
                                       AF.leadingDimension(),
                                       piv.data(),
                                       getF77Char(equed),
                                       r.data(),
                                       c.data(),
                                       B.data(),
                                       B.leadingDimension(),
                                       X.data(),
                                       X.leadingDimension(),
                                       rCond,
                                       fErr.data(),
                                       bErr.data(),
                                       work.data(),
                                       iwork.data());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
template <typename MA, typename MAF, typename VPIV, typename VR, typename VC,
          typename MB, typename MX, typename RCOND, typename FERR,
          typename BERR, typename VWORK, typename VIWORK>
typename GeMatrix<MA>::IndexType
svx(SVX::Fact           fact,
    Transpose           trans,
    GeMatrix<MA>        &A,
    GeMatrix<MAF>       &AF,
    DenseVector<VPIV>   &piv,
    SVX::Equilibration  equed,
    DenseVector<VR>     &r,
    DenseVector<VC>     &c,
    GeMatrix<MB>        &B,
    GeMatrix<MX>        &X,
    RCOND               &rCond,
    DenseVector<FERR>   &fErr,
    DenseVector<BERR>   &bErr,
    DenseVector<VWORK>  &work,
    DenseVector<VIWORK> &iwork)
{
    typedef typename GeMatrix<MA>::IndexType    IndexType;
//
//  Test the input parameters
//
#   ifndef NDEBUG
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MA>::NoView        A_org     = A;
    typename GeMatrix<MAF>::NoView       AF_org    = AF;
    typename DenseVector<VPIV>::NoView   piv_org   = piv;
    SVX::Equilibration                   equed_org = equed;
    typename DenseVector<VR>::NoView     r_org     = r;
    typename DenseVector<VC>::NoView     c_org     = c;
    typename GeMatrix<MB>::NoView        B_org     = B;
    typename GeMatrix<MX>::NoView        X_org     = X;
    RCOND                                rCond_org = rCond;
    typename DenseVector<FERR>::NoView   fErr_org  = fErr;
    typename DenseVector<BERR>::NoView   bErr_org  = bErr;
    typename DenseVector<VWORK>::NoView  work_org  = work;
    typename DenseVector<VIWORK>::NoView iwork_org = iwork;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::svx_impl(fact, trans, A, AF, piv, equed,
                                             r, c, B, X, rCond, fErr, bErr,
                                             work, iwork);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename GeMatrix<MA>::NoView        A_generic     = A;
    typename GeMatrix<MAF>::NoView       AF_generic    = AF;
    typename DenseVector<VPIV>::NoView   piv_generic   = piv;
    SVX::Equilibration                   equed_generic = equed;
    typename DenseVector<VR>::NoView     r_generic     = r;
    typename DenseVector<VC>::NoView     c_generic     = c;
    typename GeMatrix<MB>::NoView        B_generic     = B;
    typename GeMatrix<MX>::NoView        X_generic     = X;
    RCOND                                rCond_generic = rCond;
    typename DenseVector<FERR>::NoView   fErr_generic  = fErr;
    typename DenseVector<BERR>::NoView   bErr_generic  = bErr;
    typename DenseVector<VWORK>::NoView  work_generic  = work;
    typename DenseVector<VIWORK>::NoView iwork_generic = iwork;

    A       = A_org;
    AF      = AF_org;
    piv     = piv_org;
    equed   = equed_org;
    r       = r_org;
    c       = c_org;
    B       = B_org;
    X       = X_org;
    rCond   = rCond_org;
    fErr    = fErr_org;
    bErr    = bErr_org;
    work    = work_org;
    iwork   = iwork_org;

    IndexType _info = external::svx_impl(fact, trans, A, AF, piv, equed,
                                         r, c, B, X, rCond, fErr, bErr,
                                         work, iwork);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "A_org = " << A_org << std::endl;
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(AF_generic, AF, "AF_generic", "AF")) {
        std::cerr << "AF_org = " << AF_org << std::endl;
        std::cerr << "CXXLAPACK: AF_generic = " << AF_generic << std::endl;
        std::cerr << "F77LAPACK: AF = " << AF << std::endl;
        failed = true;
    }

    if (! isIdentical(piv_generic, piv, "piv_generic", "piv")) {
        std::cerr << "piv_org = " << piv_org << std::endl;
        std::cerr << "CXXLAPACK: piv_generic = " << piv_generic << std::endl;
        std::cerr << "F77LAPACK: piv = " << piv << std::endl;
        failed = true;
    }

    if (equed_generic!=equed) {
        std::cerr << "equed_org = " << equed_org << std::endl;
        std::cerr << "CXXLAPACK: equed_generic = "
                  << equed_generic << std::endl;
        std::cerr << "F77LAPACK: equed = " << equed << std::endl;
        failed = true;
    }

    if (! isIdentical(r_generic, r, "r_generic", "r")) {
        std::cerr << "r_org = " << r_org << std::endl;
        std::cerr << "CXXLAPACK: r_generic = " << r_generic << std::endl;
        std::cerr << "F77LAPACK: r = " << r << std::endl;
        failed = true;
    }

    if (! isIdentical(c_generic, c, "c_generic", "c")) {
        std::cerr << "c_org = " << c_org << std::endl;
        std::cerr << "CXXLAPACK: c_generic = " << c_generic << std::endl;
        std::cerr << "F77LAPACK: c = " << piv << std::endl;
        failed = true;
    }

    if (! isIdentical(B_generic, B, "B_generic", "B")) {
        std::cerr << "B_org = " << B_org << std::endl;
        std::cerr << "CXXLAPACK: B_generic = " << B_generic << std::endl;
        std::cerr << "F77LAPACK: B = " << B << std::endl;
        failed = true;
    }

    if (! isIdentical(X_generic, X, "X_generic", "X")) {
        std::cerr << "X_org = " << X_org << std::endl;
        std::cerr << "CXXLAPACK: X_generic = " << X_generic << std::endl;
        std::cerr << "F77LAPACK: X = " << X << std::endl;
        failed = true;
    }

    if (! isIdentical(rCond_generic, rCond, "rCond_generic", "rCond")) {
        std::cerr << "rCond_org = " << rCond_org << std::endl;
        std::cerr << "CXXLAPACK: rCond_generic = "
                  << rCond_generic << std::endl;
        std::cerr << "F77LAPACK: rCond = " << rCond << std::endl;
        failed = true;
    }

    if (! isIdentical(fErr_generic, fErr, "fErr_generic", "fErr")) {
        std::cerr << "fErr_org = " << fErr_org << std::endl;
        std::cerr << "CXXLAPACK: fErr_generic = " << fErr_generic << std::endl;
        std::cerr << "F77LAPACK: fErr = " << fErr << std::endl;
        failed = true;
    }

    if (! isIdentical(bErr_generic, bErr, "bErr_generic", "bErr")) {
        std::cerr << "bErr_org = " << bErr_org << std::endl;
        std::cerr << "CXXLAPACK: bErr_generic = " << bErr_generic << std::endl;
        std::cerr << "F77LAPACK: bErr = " << bErr << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "work_org = " << work_org << std::endl;
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (! isIdentical(iwork_generic, iwork, "iwork_generic", "iwork")) {
        std::cerr << "iwork_org = " << iwork_org << std::endl;
        std::cerr << "CXXLAPACK: iwork_generic = "
                  << iwork_generic << std::endl;
        std::cerr << "F77LAPACK: iwork = " << iwork << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, " info", "_info")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }

#   endif

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_SVX_TCC
