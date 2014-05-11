/*
 *   Copyright (c) 2014, Michael Lehn
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
       SUBROUTINE ZSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_STEQR_TCC
#define FLENS_LAPACK_IMPL_STEQR_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/lapack/typedefs.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename VD, typename VE, typename MZ, typename VWORK>
typename VD::IndexType
steqr_impl(STEQR::ComputeZ     compZ,
           DenseVector<VD>     &d,
           DenseVector<VE>     &e,
           GeMatrix<MZ>        &Z,
           DenseVector<VWORK>  &work)
{
    using cxxblas::pow;
    using std::abs;
    using std::max;
    using std::sqrt;

    typedef typename MZ::ElementType                T;
    typedef typename MZ::IndexType                  IndexType;
    typedef typename ComplexTrait<T>::PrimitiveType PT;

    const PT Zero(0), One(1), Two(2), Three(3);
    const T  CZero(0), COne(1);

    const IndexType  maxIt = 30;
    const IndexType  n     = d.length();

    const Underscore<IndexType> _;

    IndexType info = 0;

//
//  Quick return if possible
//
    if (n==0) {
        return info;
    }

    if (n==1) {
        if (compZ==STEQR::Tri) {
            Z(1,1) = One;
        }
        return info;
    }
//
//  Determine the unit roundoff and over/underflow thresholds.
//
    const PT eps      = lamch<PT>(Eps);
    const PT eps2     = pow(eps, 2);
    const PT safeMin  = lamch<PT>(SafeMin);
    const PT safeMax  = One / safeMin;
    const PT sSafeMax = sqrt(safeMax) / Three;
    const PT sSafeMin = sqrt(safeMin) / eps2;
//
//  Compute the eigenvalues and eigenvectors of the tridiagonal
//  matrix.
//
    if (compZ==STEQR::Tri) {
        Z = CZero;
        Z.diag(0) = COne;
    }

    const IndexType nMaxIt = n*maxIt;

    IndexType jTot = 0;
//
//  Determine where the matrix splits and choose QL or QR iteration
//  for each block, according to whether top or bottom diagonal
//  element is smaller.
//
    IndexType l1 = 1;
    IndexType l, lSv, lEnd, lEndSv, m, iScale;

    PT p, rt1, rt2, c, s, g, r, f, b, normA;

START:

    do {
        if (l1>n) {
            goto SORT;
        }
        if (l1>1) {
            e(l1-1) = Zero;
        }
        if (l1<=n-1) {
            for (m=l1; m<=n-1; ++m) {
                PT test = abs(e(m));
                if (test==Zero) {
                    break;
                }
                if (test<=(sqrt(abs(d(m)))*sqrt(abs(d(m+1))))*eps) {
                    e(m) = Zero;
                    break;
                }
            }
        } else {
            m = n;
        }

        l      = l1;
        lSv    = l;
        lEnd   = m;
        lEndSv = lEnd;
        l1     = m+1;

    } while (lEnd==l);
//
//  Scale submatrix in rows and columns L to LEND
//
    normA = lanst(InfinityNorm, d(_(l,lEnd)), e(_(l,max(l-1, lEnd-1))));
    iScale = 0;

    if (normA==Zero) {
        goto START;
    }
    if (normA>sSafeMax) {
        iScale = 1;
        lascl(LASCL::FullMatrix, 0, 0, normA, sSafeMax, d(_(l,lEnd)));
        lascl(LASCL::FullMatrix, 0, 0, normA, sSafeMax, e(_(l,lEnd-1)));
    } else if (normA<sSafeMin) {
        iScale = 2;
        lascl(LASCL::FullMatrix, 0, 0, normA, sSafeMin, d(_(l,lEnd)));
        lascl(LASCL::FullMatrix, 0, 0, normA, sSafeMin, e(_(l,lEnd-1)));
    }
//
//  Choose between QL and QR iteration
//
    if (abs(d(lEnd))<abs(d(l))) {
        lEnd = lSv;
        l    = lEndSv;
    }

    if (lEnd>l) {
//
//      QL Iteration
//
//      Look for small subdiagonal element.
//

QL_Start:

        while (true) {
            if (l!=lEnd) {
                for (m=l; m<=lEnd-1; ++m) {
                    PT test = pow(abs(e(m)), 2);
                    if (test<=(eps2*abs(d(m)))*abs(d(m+1))+safeMin) {
                        break;
                    }
                }
            } else {
                m = lEnd;
            }

            if (m<lEnd) {
                e(m) = Zero;
            }
            p = d(l);

            if (m==l) {
                break;
            }
//
//          If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
//          to compute its eigensystem.
//
            if (m==l+1) {
                if (compZ!=STEQR::No) {
                    laev2(d(l), e(l), d(l+1), rt1, rt2, c, s);
                    work(l)     = c;
                    work(n-1+l) = s;
                    lasr(Right, LASR::VariablePivot, LASR::Backward,
                         work(_(l,l)), work(_(n-1+l,n-1+l)),
                         Z(_(1,n),_(l,l+1)));
                } else {
                    lae2(d(l), e(l), d(l+1), rt1, rt2);
                }
                d(l)   = rt1;
                d(l+1) = rt2;
                e(l)   = Zero;
                l+=2;
                if (l<=lEnd) {
                    goto QL_Start;
                }
                goto UNDO_SCALING;
            }

            if (jTot==nMaxIt) {
                goto UNDO_SCALING;
            }
            ++jTot;
//
//          Form shift.
//
            g = (d(l+1)-p) / (Two*e(l));
            r = lapy2(g, One);
            g = d(m) - p + (e(l)/(g + sign(r,g)));

            s = One;
            c = One;
            p = Zero;
//
//          Inner loop
//
            for (IndexType i=m-1; i>=l; --i) {
                f = s*e(i);
                b = c*e(i);
                lartg(g, f, c, s, r);
                if (i!=m-1) {
                    e(i+1) = r;
                }
                g = d(i+1) - p;
                r = (d(i)-g)*s + Two*c*b;
                p = s*r;
                d(i+1) = g + p;
                g = c*r - b;
//
//              If eigenvectors are desired, then save rotations.
//
                if (compZ!=STEQR::No) {
                    work(i)     =  c;
                    work(n-1+i) = -s;
                }
            }
//
//          If eigenvectors are desired, then apply saved rotations.
//
            if (compZ!=STEQR::No) {
                lasr(Right, LASR::VariablePivot, LASR::Backward,
                     work(_(l,m-1)), work(_(n-1+l,n-2+m)),
                     Z(_(1,n),_(l,m)));
            }

            d(l) -= p;
            e(l) = g;
        }
//
//      Eigenvalue found.
//
        d(l) = p;

        ++l;
        if (l<=lEnd) {
            goto QL_Start;
        }
        goto UNDO_SCALING;

    } else {
//
//      QR Iteration
//
//      Look for small superdiagonal element.
//

QR_Start:

        while (true) {
            if (l!=lEnd) {
                for (m=l; m>=lEnd+1; --m) {
                    PT test = pow(abs(e(m-1)), 2);
                    if (test<=(eps2*abs(d(m)))*abs(d(m-1))+safeMin) {
                        break;
                    }
                }
            } else {
                m = lEnd;
            }

            if (m>lEnd) {
                e(m-1) = Zero;
            }
            p = d(l);

            if (m==l) {
                break;
            }
//
//          If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
//          to compute its eigensystem.
//
            if (m==l-1) {
                if (compZ!=STEQR::No) {
                    laev2(d(l-1), e(l-1), d(l), rt1, rt2, c, s);
                    work(m)     = c;
                    work(n-1+m) = s;
                    lasr(Right, LASR::VariablePivot, LASR::Forward,
                         work(_(m,m)), work(_(n-1+m,n-1+m)),
                         Z(_(1,n),_(l-1,l)));
                } else {
                    lae2(d(l-1), e(l-1), d(l), rt1, rt2);
                }
                d(l-1)  = rt1;
                d(l)    = rt2;
                e(l-1)  = Zero;
                l -= 2;
                if (l>=lEnd) {
                    goto QR_Start;
                }
                goto UNDO_SCALING;
            }

            if (jTot==nMaxIt) {
                goto UNDO_SCALING;
            }
            ++jTot;
//
//          Form shift.
//
            g = (d(l-1)-p) / (Two*e(l-1));
            r = lapy2(g, One);
            g = d(m) - p + (e(l-1)/(g + sign(r,g)));

            s = One;
            c = One;
            p = Zero;
//
//          Inner loop
//
            for (IndexType i=m; i<=l-1; ++i) {
                f = s*e(i);
                b = c*e(i);
                lartg(g, f, c, s, r);
                if (i!=m) {
                    e(i-1) = r;
                }
                g = d(i) - p;
                r = (d(i+1)-g)*s + Two*c*b;
                p = s*r;
                d(i) = g + p;
                g = c*r - b;
//
//              If eigenvectors are desired, then save rotations.
//
                if (compZ!=STEQR::No) {
                    work(i)     = c;
                    work(n-1+i) = s;
                }
            }
//
//          If eigenvectors are desired, then apply saved rotations.
//
            if (compZ!=STEQR::No) {
                lasr(Right, LASR::VariablePivot, LASR::Forward,
                     work(_(m,l-1)), work(_(n-1+m,n-2+l)),
                     Z(_(1,n),_(m,l)));
            }

            d(l)   -= p;
            e(l-1) = g;
        }
//
//      Eigenvalue found.
//
        d(l) = p;

        --l;
        if (l>=lEnd) {
            goto QR_Start;
        }
        goto UNDO_SCALING;

    }
//
//  Undo scaling if necessary
//
UNDO_SCALING:

    if (iScale==1) {
        lascl(LASCL::FullMatrix, 0, 0, sSafeMax, normA, d(_(lSv,lEndSv)));
        lascl(LASCL::FullMatrix, 0, 0, sSafeMax, normA, e(_(lSv,lEndSv-1)));
    } else if (iScale==2) {
        lascl(LASCL::FullMatrix, 0, 0, sSafeMin, normA, d(_(lSv,lEndSv)));
        lascl(LASCL::FullMatrix, 0, 0, sSafeMin, normA, e(_(lSv,lEndSv-1)));
    }
//
//  Check for no convergence to an eigenvalue after a total
//  of N*MAXIT iterations.
//
    if (jTot==nMaxIt) {
        for (IndexType i=1; i<=n-1; ++i) {
            if (e(i)!=Zero) {
                ++info;
            }
        }
        return info;
    }
    goto START;
//
//  Order eigenvalues and eigenvectors.
//

SORT:

    if (compZ==STEQR::No) {
//
//      Use Quick Sort
//
        lasrt(true, d);

    } else {
//
//      Use Selection Sort to minimize swaps of eigenvectors
//
        for (IndexType ii=2; ii<=n; ++ii) {
            IndexType i = ii-1;
            IndexType k = i;

            p = d(i);
            for (IndexType j=ii; j<=n; ++j) {
                if (d(j)<p) {
                    k = j;
                    p = d(j);
                }
            }
            if (k!=i) {
                d(k) = d(i);
                d(i) = p;
                blas::swap(Z(_,i), Z(_,k));
            }
        }
    }
    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename VD, typename VE, typename MZ, typename VWORK>
typename VD::IndexType
steqr_impl(STEQR::ComputeZ     compZ,
           DenseVector<VD>     &d,
           DenseVector<VE>     &e,
           GeMatrix<MZ>        &Z,
           DenseVector<VWORK>  &work)
{
    return cxxlapack::steqr(getF77Char(compZ),
                            d.length(),
                            d.data(),
                            e.data(),
                            Z.data(),
                            Z.leadingDimension(),
                            work.data());
}

} // namespace external

#endif

//== public interface ==========================================================

template <typename VD, typename VE, typename MZ, typename VWORK>
typename RestrictTo<IsRealDenseVector<VD>::value
                 && IsRealDenseVector<VE>::value
                 && IsComplexGeMatrix<MZ>::value
                 && IsRealDenseVector<VWORK>::value,
         typename RemoveRef<VD>::Type::IndexType>::Type
steqr(STEQR::ComputeZ  compZ,
      VD               &&d,
      VE               &&e,
      MZ               &&Z,
      VWORK            &&work)
{
    LAPACK_DEBUG_OUT("steqr");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<VD>::Type     VectorD;
    typedef typename VectorD::IndexType      IndexType;

#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<VE>::Type     VectorE;
    typedef typename RemoveRef<MZ>::Type     MatrixZ;
    typedef typename RemoveRef<VWORK>::Type  VectorWork;
#   endif

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(d.firstIndex()==1);
    ASSERT(e.firstIndex()==1);
    ASSERT(Z.firstRow()==1);
    ASSERT(Z.firstCol()==1);

    const IndexType  n = d.length();
    ASSERT(e.length()==n-1);

    if (compZ==STEQR::Orig ||compZ==STEQR::Tri) {
        ASSERT(Z.numRows()==n);
        ASSERT(Z.numCols()==n);
    }

    ASSERT(work.firstIndex()==1);
    ASSERT(work.length()==std::max(1,2*n-2));
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename VectorD::NoView      d_org    = d;
    typename VectorE::NoView      e_org    = e;
    typename MatrixZ::NoView      Z_org    = Z;
    typename VectorWork::NoView   work_org = work;
#   endif

//
//  Call implementation
//
    const IndexType info = LAPACK_SELECT::steqr_impl(compZ, d, e, Z, work);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename VectorD::NoView      d_generic    = d;
    typename VectorE::NoView      e_generic    = e;
    typename MatrixZ::NoView      Z_generic    = Z;
    typename VectorWork::NoView   work_generic = work;

    d    = d_org;
    e    = e_org;
    Z    = Z_org;
    work = work_org;

    const IndexType _info = external::steqr_impl(compZ, d, e, Z, work);

    bool failed = false;
    if (! isIdentical(_info, info, "_info", "info")) {
        std::cerr << "CXXLAPACK: _info = " << _info << std::endl;
        std::cerr << "F77LAPACK: info = " << info << std::endl;
        failed = true;
    }
    if (! isIdentical(d_generic, d, "d_generic", "d")) {
        std::cerr << "CXXLAPACK: d_generic = " << d_generic << std::endl;
        std::cerr << "F77LAPACK: d = " << d << std::endl;
        failed = true;
    }
    if (! isIdentical(e_generic, e, "e_generic", "e")) {
        std::cerr << "CXXLAPACK: e_generic = " << e_generic << std::endl;
        std::cerr << "F77LAPACK: e = " << e << std::endl;
        failed = true;
    }
    if (! isIdentical(Z_generic, Z, "Z_generic", "Z")) {
        std::cerr << "CXXLAPACK: Z_generic = " << Z_generic << std::endl;
        std::cerr << "F77LAPACK: Z = " << Z << std::endl;
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

    return info;
}


} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_STEQR_TCC
