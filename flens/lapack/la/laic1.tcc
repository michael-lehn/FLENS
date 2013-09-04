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
       SUBROUTINE DLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C )
 *
 *  -- LAPACK auxiliary routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 *
 */

#ifndef FLENS_LAPACK_LA_LAIC1_TCC
#define FLENS_LAPACK_LA_LAIC1_TCC 1

#include <flens/lapack/lapack.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- laic1 [real variant] ------------------------------------------------------

template <typename VX, typename SEST, typename VW, typename GAMMA,
          typename SESTPR, typename S, typename C>
void
laic1_impl(LAIC1::Job               job,
           const DenseVector<VX>    &x,
           SEST                     sEst,
           const DenseVector<VW>    &w,
           const GAMMA              &gamma,
           SESTPR                   &sEstPr,
           S                        &s,
           C                        &c)
{
    using std::abs;
    using std::max;
    using std::sqrt;
    using cxxblas::conjugate;
    using cxxblas::real;

    typedef typename DenseVector<VX>::ElementType              ElementType;
    typedef typename ComplexTrait<ElementType>::PrimitiveType  PrimitiveType;

    const PrimitiveType Zero(0), One(1), Two(2), Half(0.5), Four(4);

    const PrimitiveType eps   = lamch<PrimitiveType>(Eps);
    const ElementType   alpha = blas::dot(x, w);

    const PrimitiveType absAlpha = abs(alpha);
    const PrimitiveType absGamma = abs(gamma);
    const PrimitiveType absEst   = abs(sEst);

    PrimitiveType s1, s2, tmp, zeta1, zeta2, b, t, scl, normA;
    ElementType   cosine, sine;

    if (job==LAIC1::Max) {
//
//      Estimating largest singular value
//
//      special cases
//
        if (sEst==Zero) {
            s1 =max(absGamma, absAlpha);
            if (s1==Zero) {
                s = Zero;
                c = One;
                sEstPr = Zero;
            } else {
                s = alpha / s1;
                c = gamma / s1;
                tmp = sqrt(cxxblas::real(s*conjugate(s)) + cxxblas::real(c*conjugate(c)));
                s /= tmp;
                c /= tmp;
                sEstPr = s1*tmp;
            }
            return;
        } else if (absGamma<=eps*absEst) {
            s = One;
            c = Zero;
            tmp = max(absEst, absAlpha);
            s1 = absEst / tmp;
            s2 = absAlpha / tmp;
            sEstPr = tmp*sqrt(s1*s1 + s2*s2);
            return;
        } else if (absAlpha<=eps*absEst) {
            s1 = absGamma;
            s2 = absEst;
            if (s1<=s2) {
                s = One;
                c = Zero;
                sEstPr = s2;
            } else {
                s = Zero;
                c = One;
                sEstPr = s1;
            }
            return;
        } else if (absEst<=eps*absAlpha || absEst<=eps*absGamma) {
            s1 = absGamma;
            s2 = absAlpha;
            if (s1<=s2) {
                tmp = s1 / s2;
                if (IsComplex<ElementType>::value) {
                    scl = sqrt(One + tmp*tmp);
                    sEstPr = s2*scl;
                    s = (alpha/s2) / scl;
                    c = (gamma/s2) / scl;
                } else {
                    s = sqrt(One + tmp*tmp);
                    sEstPr = cxxblas::real(s2*s);
                    c = (gamma/s2) / s;
                    s = sign(One, cxxblas::real(alpha)) / s;
                }
            } else {
                tmp = s2 / s1;
                if (IsComplex<ElementType>::value) {
                    scl = sqrt(One + tmp*tmp);
                    sEstPr = s1*scl;
                    s = (alpha/s1) / scl;
                    c = (gamma/s1) / scl;
                } else {
                    c = sqrt(One + tmp*tmp);
                    sEstPr = cxxblas::real(s1*c);
                    s = (alpha/s1) / c;
                    c = sign(One, cxxblas::real(gamma)) / c;
                }
            }
            return;
        } else {
//
//          normal case
//
            if (IsComplex<ElementType>::value) {
                zeta1 = absAlpha / absEst;
                zeta2 = absGamma / absEst;
            } else {
                zeta1 = cxxblas::real(alpha / absEst);
                zeta2 = cxxblas::real(gamma / absEst);
            }

            b = (One - zeta1*zeta1 - zeta2*zeta2)*Half;
            c = zeta1*zeta1;
            if (b>Zero) {
                t = cxxblas::real(c) / (b+sqrt(b*b+cxxblas::real(c)));
            } else {
                t = sqrt(b*b+cxxblas::real(c)) - b;
            }

            if (IsComplex<ElementType>::value) {
                sine   = -(alpha / absEst) / t;
                cosine = -(gamma / absEst) / (One +t);
            } else {
                sine   = -zeta1 / t;
                cosine = -zeta2 / (One +t);
            }
            tmp = sqrt(cxxblas::real(sine*conjugate(sine)) + cxxblas::real(cosine*conjugate(cosine)));
            s = sine / tmp;
            c = cosine / tmp;
            sEstPr = sqrt(t + One)*absEst;
            return;
        }

    } else if (job==LAIC1::Min) {
//
//      Estimating smallest singular value
//
//      special cases
//
        if (sEst==Zero) {
            sEstPr = Zero;
            if (max(absGamma, absAlpha)==Zero) {
                sine   = One;
                cosine = Zero;
            } else {
                sine   = -gamma;
                cosine = alpha;
            }
            s1 = max(abs(sine), abs(cosine));
            s = sine / s1;
            c = cosine / s1;
            tmp = sqrt(cxxblas::real(s*conjugate(s)) + cxxblas::real(c*conjugate(c)));
            s /= tmp;
            c /= tmp;
            return;
        } else if (absGamma<=eps*absEst) {
            s = Zero;
            c = One;
            sEstPr = absGamma;
            return;
        } else if (absAlpha<=eps*absEst) {
            s1 = absGamma;
            s2 = absEst;
            if (s1<=s2) {
                s = Zero;
                c = One;
                sEstPr = s1;
            } else {
                s = One;
                c = Zero;
                sEstPr = s2;
            }
            return;
        } else if (absEst<=eps*absAlpha || absEst<=eps*absGamma) {
            s1 = absGamma;
            s2 = absAlpha;
            if (s1<=s2) {
                tmp = s1 / s2;
                if (IsComplex<ElementType>::value) {
                    scl = sqrt(One + tmp*tmp);
                    sEstPr = absEst*(tmp/scl);
                    s = -(conjugate(gamma) / s2) / scl;
                    c =  (conjugate(alpha) / s2) / scl;
                } else {
                    c = sqrt(One + tmp*tmp);
                    sEstPr = absEst*(tmp/cxxblas::real(c));
                    s = -(gamma / s2) / c;
                    c = sign(One, cxxblas::real(alpha)) / c;
                }
            } else {
                tmp = s2 / s1;
                if (IsComplex<ElementType>::value) {
                    scl = sqrt(One + tmp*tmp);
                    sEstPr = absEst / scl;
                    s = -(conjugate(gamma) / s1) / scl;
                    c =  (conjugate(alpha) / s1) / scl;
                } else {
                    s = sqrt(One + tmp*tmp);
                    sEstPr = absEst / cxxblas::real(s);
                    c = (alpha / s1) / s;
                    s = -sign(One, cxxblas::real(gamma)) / s;
                }
            }
            return;
        } else {
//
//          normal case
//
            if (IsComplex<ElementType>::value) {
                zeta1 = absAlpha / absEst;
                zeta2 = absGamma / absEst;
            } else {
                zeta1 = cxxblas::real(alpha) / absEst;
                zeta2 = cxxblas::real(gamma) / absEst;
            }

            if (IsComplex<ElementType>::value) {
                normA = max(One + zeta1*zeta1 + zeta1*zeta2,
                            zeta1*zeta2 + zeta2*zeta2);
            } else {
                normA = max(One + zeta1*zeta1 + abs(zeta1*zeta2),
                            abs(zeta1*zeta2) + zeta2*zeta2);
            }
//
//          See if root is closer to zero or to ONE
//
            PrimitiveType test = One + Two*(zeta1-zeta2)*(zeta1+zeta2);
            if (test>=Zero) {
//
//              root is close to zero, compute directly
//
                b = (zeta1*zeta1 + zeta2*zeta2 + One)*Half;
                c = zeta2*zeta2;
                t = cxxblas::real(c) / (b + sqrt(abs(b*b-cxxblas::real(c))));
                if (IsComplex<ElementType>::value) {
                    sine   = (alpha/absEst) / (One - t);
                    cosine = -(gamma/absEst) / t;
                } else {
                    sine   = zeta1 / (One - t);
                    cosine = -zeta2 / t;
                }
                sEstPr = sqrt(t + Four*eps*eps*normA) * absEst;
            } else {
//
//              root is closer to ONE, shift by that amount
//
                b = (zeta2*zeta2 + zeta1*zeta1 - One)*Half;
                c = zeta1*zeta1;
                if (b>=Zero) {
                    t = -cxxblas::real(c) / (b + sqrt(b*b + cxxblas::real(c)));
                } else {
                    t = b - sqrt(b*b + cxxblas::real(c));
                }
                if (IsComplex<ElementType>::value) {
                    sine   = -(alpha/absEst) / t;
                    cosine = -(gamma/absEst) / (One + t);
                } else {
                    sine   = -zeta1 / t;
                    cosine = -zeta2 / (One + t);
                }
                sEstPr = sqrt(One + t + Four*eps*eps*normA)*absEst;
            }
            tmp = sqrt(cxxblas::real(sine*conjugate(sine)) + cxxblas::real(cosine*conjugate(cosine)));
            s = sine / tmp;
            c = cosine / tmp;
            return;

        }
    }

}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- laic1 [real variant] ------------------------------------------------------

template <typename VX, typename SEST, typename VW, typename GAMMA,
          typename SESTPR, typename S, typename C>
void
laic1_impl(LAIC1::Job               job,
           const DenseVector<VX>    &x,
           SEST                     sEst,
           const DenseVector<VW>    &w,
           const GAMMA              &gamma,
           SESTPR                   &sEstPr,
           S                        &s,
           C                        &c)
{
    typedef typename DenseVector<VX>::IndexType  IndexType;

    cxxlapack::laic1<IndexType>(IndexType(job),
                                x.length(),
                                x.data(),
                                sEst,
                                w.data(),
                                gamma,
                                sEstPr,
                                s,
                                c);
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- laic1 [real variant] ------------------------------------------------------

template <typename VX, typename SEST, typename VW, typename GAMMA,
          typename SESTPR, typename S, typename C>
void
laic1(LAIC1::Job                job,
      const DenseVector<VX>     &x,
      SEST                      sEst,
      const DenseVector<VW>     &w,
      const GAMMA               &gamma,
      SESTPR                    &sEstPr,
      S                         &s,
      C                         &c)
{

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(x.length()==w.length());
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    SESTPR   sEstPr_org = sEstPr;
    S        s_org      = s;
    C        c_org      = c;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::laic1_impl(job, x, sEst, w, gamma, sEstPr, s, c);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    SESTPR   sEstPr_generic = sEstPr;
    S        s_generic      = s;
    C        c_generic      = c;

//
//  restore output arguments
//
    sEstPr = sEstPr_org;
    s      = s_org;
    c      = c_org;

//
//  Compare generic results with results from the native implementation
//
    external::laic1_impl(job, x, sEst, w, gamma, sEstPr, s, c);

    bool failed = false;
    if (! isIdentical(sEstPr_generic, sEstPr, "sEstPr_generic", "sEstPr")) {
        std::cerr << "CXXLAPACK: sEstPr_generic = "
                  << sEstPr_generic << std::endl;
        std::cerr << "F77LAPACK: sEstPr = " << sEstPr << std::endl;
        failed = true;
    }
    if (! isIdentical(s_generic, s, "s_generic", "s")) {
        std::cerr << "CXXLAPACK: s_generic = " << s_generic << std::endl;
        std::cerr << "F77LAPACK: s = " << s << std::endl;
        failed = true;
    }
    if (! isIdentical(c_generic, c, "c_generic", "c")) {
        std::cerr << "CXXLAPACK: c_generic = " << c_generic << std::endl;
        std::cerr << "F77LAPACK: c = " << c << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: laic1.tcc" << std::endl;
        //ASSERT(0);
    } else {
        // std::cerr << "passed: laic1.tcc" << std::endl;
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAIC1_TCC
