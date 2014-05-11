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
       SUBROUTINE DSTERF( N, D, E, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_STERF_TCC
#define FLENS_LAPACK_IMPL_STERF_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/lapack/typedefs.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename VD, typename VE>
typename VD::IndexType
sterf_impl(DenseVector<VD> &d,
           DenseVector<VE> &e)
{
    using cxxblas::pow;
    using std::abs;
    using std::sqrt;

    typedef typename VD::ElementType  T;
    typedef typename VD::IndexType    IndexType;

    const T Zero(0), One(1), Two(2), Three(3);

    const Underscore<IndexType>   _;

    const IndexType n     = d.length();
    const IndexType maxIt = 30;

    IndexType info = 0;

//
//  Quick return if possible
//
    if (n<=1) {
        return info;
    }
//
//  Determine the unit roundoff for this environment.
//
    const T eps      = lamch<T>(Precision);
    const T eps2     = pow(eps, 2);
    const T safeMin  = lamch<T>(SafeMin);
    const T safeMax  = One / SafeMin;
    const T sSaveMax = sqrt(safeMax) / Three;
    const T sSaveMin = sqrt(safeMin) / eps2;
    //const T rMax     = lamch<T>(OverflowThreshold);
//
//  Compute the eigenvalues of the tridiagonal matrix.
//
    const IndexType nMaxIt = n*maxIt;

    T   sigma = Zero;

    IndexType jTot  = 0;
//
//  Determine where the matrix splits and choose QL or QR iteration
//  for each block, according to whether top or bottom diagonal
//  element is smaller.
//
    IndexType l1 = 1;
    IndexType l, lSv, lEnd, lEndSv, m, iScale;

    T normA, p, rte, rt1, rt2, r, c, s, alpha, gamma, bb, oldC, oldGamma;

START:

    do {
        if (l1>n) {
            goto SORT;
        }
        if (l1>1) {
            e(l1-1) = Zero;
        }
        for (m=l1; m<=n-1; ++m) {
            if (abs(e(m))<=(sqrt(abs(d(m)))*sqrt(abs(d(m+1))))*eps) {
                e(m) = Zero;
                break;
            }
        }

        l      = l1;
        lSv    = l;
        lEnd   = m;
        lEndSv = lEnd;
        l1     = m + 1;

    } while (lEnd==l);

//
//  Scale submatrix in rows and columns L to LEND
//
    normA = lanst(MaximumNorm, d(_(l,lEnd)), e(_(l,lEnd-1)));
    iScale = 0;
    if (normA==Zero) {
        goto START;
    }
    if (normA>sSaveMax) {
        iScale = 1;

        lascl(LASCL::FullMatrix, 0, 0, normA, sSaveMax, d(_(l,lEnd)));
        lascl(LASCL::FullMatrix, 0, 0, normA, sSaveMax, e(_(l,lEnd-1)));
    } else if (normA<sSaveMin) {
        iScale = 2;
        lascl(LASCL::FullMatrix, 0, 0, normA, sSaveMin, d(_(l,lEnd)));
        lascl(LASCL::FullMatrix, 0, 0, normA, sSaveMin, e(_(l,lEnd-1)));
    }

    for (IndexType i=l; i<=lEnd-1; ++i) {
        e(i) = pow(e(i), 2);
    }
//
//  Choose between QL and QR iteration
//
    if (abs(d(lEnd))<abs(d(l))) {
        lEnd = lSv;
        l    = lEndSv;
    }

    if (lEnd>=l) {
//
//      QL Iteration
//
//      Look for small subdiagonal element.
//
QL_Start:

        while (true) {
            if (l!=lEnd) {
                for (IndexType m=l; l<=lEnd-1; ++l) {
                    if (abs(e(m))<=eps2*abs(d(m)*d(m+1))) {
                        break;
                    }
                }
            }

            if (m<lEnd) {
                e(m) = Zero;
            }
            p = d(l);
            if (m==l) {
                break;
            }
//
//          If remaining matrix is 2 by 2, use DLAE2 to compute its
//          eigenvalues.
//
            if (m==l+1) {
                rte = sqrt(e(l));

                lae2(d(l), rte, d(l+1), rt1, rt2);

                d(l)   = rt1;
                d(l+1) = rt2;
                e(l)   = Zero;
                l += 2;
                if (l<=lEnd) {
                    continue;
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
            T rte = sqrt(e(l));
            sigma = (d(l+1)-p) / (Two*rte);
            r     = lapy2(sigma, One);
            sigma = p - (rte / (sigma+sign(r,sigma)));

            c     = One;
            s     = Zero;
            gamma = d(m) - sigma;
            p     = gamma*gamma;
//
//          Inner loop
//
            for (IndexType i=m-1; i>=l; --i) {
                bb = e(i);
                r  = p + bb;
                if (i!=m-1) {
                    e(i+1) = s*r;
                }
                oldC = c;
                c = p / r;
                s = bb / r;
                oldGamma = gamma;
                alpha = d(i);
                gamma = c*(alpha-sigma) - s*oldGamma;
                d(i+1) = oldGamma + (alpha-gamma);
                if (c!=Zero) {
                    p = (gamma*gamma)/c;
                } else {
                    p = oldC*bb;
                }
            }

            e(l) = s*p;
            d(l) = sigma + gamma;
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
            for (m=l; m>=lEnd+1; --m) {
                if (abs(e(m-1))<=eps2*abs(d(m)*d(m-1))) {
                    break;
                }
            }

            if (m>lEnd) {
                e(m-1) = Zero;
            }
            p = d(l);
            if (m==l) {
                break;
            }
//
//          If remaining matrix is 2 by 2, use DLAE2 to compute its
//          eigenvalues.
//
            if (m==l-1) {
                rte = sqrt(e(l-1));

                lae2(d(l), rte, d(l-1), rt1, rt2);

                d(l)   = rt1;
                d(l-1) = rt2;
                e(l-1) = Zero;

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
            rte   = sqrt(e(l-1));
            sigma = (d(l-1)-p) / (Two*rte);
            r     = lapy2(sigma, One);
            sigma = p - (rte / (sigma+sign(r, sigma)));

            c     = One;
            s     = Zero;
            gamma = d(m) - sigma;
            p     = gamma*gamma;
//
//          Inner loop
//
            for (IndexType i=m; i<=l-1; ++i) {
                bb = e(i);
                r  = p + bb;
                if (i!=m) {
                    e(i-1) = s*r;
                }
                oldC = c;
                c = p / r;
                s = bb / r;
                oldGamma = gamma;
                alpha = d(i+1);
                gamma = c*(alpha-sigma) - s*oldGamma;
                d(i) = oldGamma + (alpha-gamma);
                if (c!=Zero) {
                    p = (gamma*gamma) / c;
                } else {
                    p = oldC*bb;
                }
            }

            e(l-1) = s*p;
            d(l)   = sigma + gamma;
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
        lascl(LASCL::FullMatrix, 0, 0, sSaveMax, normA, d(_(lSv,lEndSv)));
    }
    if (iScale==2) {
        lascl(LASCL::FullMatrix, 0, 0, sSaveMin, normA, d(_(lSv,lEndSv)));
    }
//
//  Check for no convergence to an eigenvalue after a total
//  of N*MAXIT iterations.
//
    if (jTot<nMaxIt) {
        goto START;
    }
    for (IndexType i=1; i<=n-1; ++i) {
        if (e(i)!=Zero) {
            ++info;
        }
    }
    goto RETURN;
//
//  Sort eigenvalues in increasing order.
//
SORT:

    lasrt(true, d);

RETURN:

    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename VD, typename VE>
typename VD::IndexType
sterf_impl(DenseVector<VD> &d,
           DenseVector<VE> &e)
{
    return cxxlapack::sterf(d.length(),
                            d.data(),
                            e.data());
}

} // namespace external

#endif

//== public interface ==========================================================

template <typename VD, typename VE>
typename RestrictTo<IsRealDenseVector<VD>::value
                 && IsRealDenseVector<VE>::value,
         typename RemoveRef<VD>::Type::IndexType>::Type
sterf(VD  &&d,
      VE  &&e)
{
    LAPACK_DEBUG_OUT("sterf");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<VD>::Type    VectorD;
    typedef typename VectorD::IndexType     IndexType;

#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<VE>::Type    VectorE;
#   endif

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(d.firstIndex()==1);
    ASSERT(e.firstIndex()==1);

    const IndexType  n = d.length();
    ASSERT(e.length()==n-1);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename VectorD::NoView      d_org    = d;
    typename VectorE::NoView      e_org    = e;
#   endif

//
//  Call implementation
//
    const IndexType info = LAPACK_SELECT::sterf_impl(d, e);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename VectorD::NoView      d_generic    = d;
    typename VectorE::NoView      e_generic    = e;

    d = d_org;
    e = e_org;

    const IndexType _info = external::sterf_impl(d, e);

    bool failed = false;
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
    if (! isIdentical(info, _info, "info", "_info")) {
        std::cerr << "CXXLAPACK: info= " << info << std::endl;
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

#endif // FLENS_LAPACK_IMPL_STERF_TCC
