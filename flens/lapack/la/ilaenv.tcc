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
 *    INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
 *
 *    http://www.netlib.org/lapack/util/ilaenv.f
 *
 *  -- LAPACK auxiliary routine (version 3.2.1)                        --
 *
 *  -- April 2009                                                      --
 *
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *
 */

#ifndef FLENS_LAPACK_LA_ILAENV_TCC
#define FLENS_LAPACK_LA_ILAENV_TCC 1

#include <cxxstd/complex.h>
#include <cxxstd/cstring.h>

#include <flens/auxiliary/auxiliary.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <int n>
bool
isSame(const char *a, const char *b)
{
    for (int i=0; i<n; ++i) {
        if (a[i]!=b[i]) {
            return false;
        }
    }
    return true;
}

namespace generic {

template <typename T>
int
ilaenv_impl(int spec, const char *name, const char *opts,
            int n1, int n2, int n3, int n4)
{
    int result = -1;
    const char *c2 = name;
    const char *c3 = name + 2;
    const char *c4 = name + 3;

    int nb, nbMin, nx;

    switch (spec) {
//
//      optimal blocksize
//
        case 1:
            nb = 32;
            if (isSame<2>(c2, "GE")) {
                if (isSame<3>(c3, "TRF")) {
                    if (IsNotComplex<T>::value) {
                        nb = 64;
                    } else {
                        nb = 64;
                    }
                } else if (isSame<3>(c3, "QRF") || isSame<3>(c3, "RQF")
                       ||  isSame<3>(c3, "LQF") || isSame<3>(c3, "QLF")) {
                    if (IsNotComplex<T>::value) {
                        nb = 32;
                    } else {
                        nb = 32;
                    }
                } else if (isSame<3>(c3, "HRD")) {
                    if (IsNotComplex<T>::value) {
                        nb = 32;
                    } else {
                        nb = 32;
                    }
                } else if (isSame<3>(c3, "BRD")) {
                    if (IsNotComplex<T>::value) {
                        nb = 32;
                    } else {
                        nb = 32;
                    }
                } else if (isSame<3>(c3, "TRI")) {
                    if (IsNotComplex<T>::value) {
                        nb = 64;
                    } else {
                        nb = 64;
                    }
                }
            } else if (isSame<2>(c2, "PO")) {
                if (isSame<3>(c3, "TRF")) {
                    if (IsNotComplex<T>::value) {
                        nb = 64;
                    } else {
                        nb = 64;
                    }
                }
            } else if (isSame<2>(c2, "SY")) {
                if (isSame<3>(c3, "TRF")) {
                    if (IsNotComplex<T>::value) {
                        nb = 64;
                    } else {
                        nb = 64;
                    }
                } else if (IsNotComplex<T>::value && isSame<3>(c3, "TRD")) {
                    nb = 32;
                } else if (IsNotComplex<T>::value && isSame<3>(c3, "GST")) {
                    nb = 64;
                }
            } else if (IsComplex<T>::value &&  isSame<2>(c2, "HE")) {
                if (isSame<3>(c3, "TRF")) {
                    nb = 64;
                } else if (isSame<3>(c3, "TRD")) {
                    nb = 32;
                } else if (isSame<3>(c3, "GST")) {
                    nb = 64;
                }
            } else if (IsNotComplex<T>::value && isSame<2>(c2, "OR")) {
                if (isSame<1>(c3, "G")) {
                    if (isSame<2>(c4, "QR") || isSame<2>(c4, "RQ")
                     || isSame<2>(c4, "LQ") || isSame<2>(c4, "QL")
                     || isSame<2>(c4, "HR") || isSame<2>(c4, "TR")
                     || isSame<2>(c4, "BR"))
                    {
                        nb = 32;
                    }
                } else if (isSame<1>(c3, "M")) {
                    if (isSame<2>(c4, "QR") || isSame<2>(c4, "RQ")
                     || isSame<2>(c4, "LQ") || isSame<2>(c4, "QL")
                     || isSame<2>(c4, "HR") || isSame<2>(c4, "TR")
                     || isSame<2>(c4, "BR"))
                    {
                        nb = 32;
                    }
                }
            } else if (IsComplex<T>::value && isSame<2>(c2, "UN")) {
                if (isSame<1>(c3, "G")) {
                    if (isSame<2>(c4, "QR") || isSame<2>(c4, "RQ")
                     || isSame<2>(c4, "LQ") || isSame<2>(c4, "QL")
                     || isSame<2>(c4, "HR") || isSame<2>(c4, "TR")
                     || isSame<2>(c4, "BR"))
                    {
                        nb = 32;
                    }
               } else if (isSame<1>(c3, "M")) {
                    if (isSame<2>(c4, "QR") || isSame<2>(c4, "RQ")
                     || isSame<2>(c4, "LQ") || isSame<2>(c4, "QL")
                     || isSame<2>(c4, "HR") || isSame<2>(c4, "TR")
                     || isSame<2>(c4, "BR"))
                    {
                        nb = 32;
                    }
               }
            } else if (isSame<2>(c2, "GB")) {
               if (isSame<3>(c3, "TRF")) {
                  if (IsNotComplex<T>::value) {
                     if (n4<=64) {
                        nb = 1;
                     } else {
                        nb = 32;
                     }
                  } else {
                     if (n4<=64) {
                        nb = 1;
                     } else {
                        nb = 32;
                     }
                  }
               }
            } else if (isSame<2>(c2, "PB")) {
               if (isSame<3>(c3, "TRF")) {
                  if (IsNotComplex<T>::value) {
                     if (n2<=64) {
                        nb = 1;
                     } else {
                        nb = 32;
                     }
                  } else {
                     if (n2<=64 ) {
                        nb = 1;
                     } else {
                        nb = 32;
                     }
                  }
               }
            } else if (isSame<2>(c2, "TR")) {
               if (isSame<3>(c3, "TRI")) {
                  if (IsNotComplex<T>::value) {
                     nb = 64;
                  } else {
                     nb = 64;
                  }
               }
            } else if (isSame<2>(c2, "LA")) {
               if (isSame<3>(c3, "UUM")) {
                  if (IsNotComplex<T>::value) {
                     nb = 64;
                  } else {
                     nb = 64;
                  }
               }
            } else if (IsNotComplex<T>::value && isSame<2>(c2, "ST") ) {
               if (isSame<3>(c3, "EBZ")) {
                  nb = 1;
               }
            }
            result = nb;
            break;
//
//      minimal blocksize
//
        case 2:
            nbMin = 2;
            if (isSame<2>(c2, "GE")) {
                if (isSame<3>(c3, "QRF") || isSame<3>(c3, "RQF")
                 || isSame<3>(c3, "LQF") || isSame<3>(c3, "QLF")) {
                    if (IsNotComplex<T>::value ) {
                         nbMin = 2;
                    } else {
                        nbMin = 2;
                    }
                } else if (isSame<3>(c3, "HRD")) {
                    if (IsNotComplex<T>::value ) {
                        nbMin = 2;
                    } else {
                        nbMin = 2;
                    }
                } else if (isSame<3>(c3, "BRD")) {
                    if (IsNotComplex<T>::value ) {
                        nbMin = 2;
                    } else {
                        nbMin = 2;
                    }
                } else if (isSame<3>(c3, "TRI")) {
                    if (IsNotComplex<T>::value ) {
                        nbMin = 2;
                    } else {
                        nbMin = 2;
                    }
                }
            } else if (isSame<2>(c2, "SY")) {
                if (isSame<3>(c3, "TRF")) {
                    if (IsNotComplex<T>::value ) {
                        nbMin = 8;
                    } else {
                        nbMin = 8;
                    }
                } else if (IsNotComplex<T>::value && isSame<3>(c3, "TRD")) {
                    nbMin = 2;
                }
            } else if (IsComplex<T>::value && isSame<2>(c2, "HE")) {
                if (isSame<3>(c3, "TRD")) {
                    nbMin = 2;
                }
            } else if (IsNotComplex<T>::value && isSame<2>(c2, "OR")) {
                if (isSame<1>(c3, "G")) {
                    if (isSame<2>(c4, "QR") || isSame<2>(c4, "RQ")
                     || isSame<2>(c4, "LQ") || isSame<2>(c4, "QL")
                     || isSame<2>(c4, "HR") || isSame<2>(c4, "TR")
                     || isSame<2>(c4, "BR")) {
                        nbMin = 2;
                    }
                } else if (isSame<1>(c3, "M")) {
                    if (isSame<2>(c4, "QR") || isSame<2>(c4, "RQ")
                     || isSame<2>(c4, "LQ") || isSame<2>(c4, "QL")
                     || isSame<2>(c4, "HR") || isSame<2>(c4, "TR")
                     || isSame<2>(c4, "BR")) {
                        nbMin = 2;
                    }
                }
            } else if (IsComplex<T>::value && isSame<2>(c2, "UN")) {
                if (isSame<1>(c3, "G")) {
                    if (isSame<2>(c4, "QR") || isSame<2>(c4, "RQ")
                     || isSame<2>(c4, "LQ") || isSame<2>(c4, "QL")
                     || isSame<2>(c4, "HR") || isSame<2>(c4, "TR")
                     || isSame<2>(c4, "BR")) {
                        nbMin = 2;
                    }
                } else if (isSame<1>(c3, "M")) {
                    if (isSame<2>(c4, "QR") || isSame<2>(c4, "RQ")
                     || isSame<2>(c4, "LQ") || isSame<2>(c4, "QL")
                     || isSame<2>(c4, "HR") || isSame<2>(c4, "TR")
                     || isSame<2>(c4, "BR")) {
                        nbMin = 2;
                    }
                }
            }
            result = nbMin;
            break;
//
//      crossover point
//
        case 3:
            nx = 32;
            if (isSame<2>(c2, "GE")) {
                if (isSame<3>(c3, "QRF") || isSame<3>(c3, "RQF")
                 || isSame<3>(c3, "LQF") || isSame<3>(c3, "QLF")) {
                    if (IsNotComplex<T>::value) {
                        nx =  128;
                    } else {
                        nx =  128;
                    }
                } else if (isSame<3>(c3, "HRD")) {
                    if (IsNotComplex<T>::value) {
                        nx =  128;
                    } else {
                        nx =  128;
                    }
                } else if (isSame<3>(c3, "BRD")) {
                    if (IsNotComplex<T>::value) {
                        nx =  128;
                    } else {
                        nx =  128;
                    }
                }
            } else if (isSame<2>(c2, "SY")) {
                if (IsNotComplex<T>::value && isSame<3>(c3, "TRD")) {
                    nx =  32;
                }
            } else if (IsComplex<T>::value && isSame<2>(c2, "HE")) {
                if (isSame<3>(c3, "TRD")) {
                    nx =  32;
                }
            } else if (IsNotComplex<T>::value && isSame<2>(c2, "OR")) {
                if (isSame<1>(c3, "G")) {
                    if (isSame<2>(c4, "QR") || isSame<2>(c4, "RQ")
                     || isSame<2>(c4, "LQ") || isSame<2>(c4, "QL")
                     || isSame<2>(c4, "HR") || isSame<2>(c4, "TR")
                     || isSame<2>(c4, "BR")) {
                        nx =  128;
                    }
                }
            } else if (IsComplex<T>::value && isSame<2>(c2, "UN")) {
                if (isSame<1>(c3, "G")) {
                    if (isSame<2>(c4, "QR") || isSame<2>(c4, "RQ")
                     || isSame<2>(c4, "LQ") || isSame<2>(c4, "QL")
                     || isSame<2>(c4, "HR") || isSame<2>(c4, "TR")
                     || isSame<2>(c4, "BR")) {
                        nx =  128;
                    }
                }
            }
            result = static_cast<int>(static_cast<double>(nx)/1.5);
            break;
//
//      DEPRECATED: number of shifts used in the nonsymmetric eigenvalue
//                  routines
//
        case 4:
            ASSERT(0);
            break;
//
//      12 <= pec<= 16: hseqr or one of its subroutines ..
//
        case 12:
        case 13:
        case 14:
        case 15:
        case 16:
            result = iparmq<T>(spec, name, opts, n1, n2, n3, n4);
            break;

        default:
            ASSERT(0);
            result = -1;
            break;

    }

#   ifdef LOG_ILAENV
    std::cerr << "ILAENV_GENERIC: "
              << "spec = " << spec
              << ", name = " << name
              << ", opts = " << opts
              << ", n1 " << n1
              << ", n2 " << n2
              << ", n3 " << n3
              << ", n4 " << n4
              << ", result = " << result
              << std::endl;
#   endif


    return result;
}

} // namespace generic

//== interface for native lapack ===============================================

#if defined CHECK_CXXLAPACK || defined USE_NATIVE_ILAENV

#ifndef LAPACK_DECL
#   define  LAPACK_DECL(x)    x##_

#endif

extern "C" {

INTEGER
LAPACK_DECL(ilaenv)(const INTEGER *SPEC,
                    const char *NAME,
                    const char *OPTS,
                    const INTEGER *N1,
                    const INTEGER *N2,
                    const INTEGER *N3,
                    const INTEGER *N4,
                    int NAME_LEN,
                    int OPTS_LEN);

} // extern "C"



template <typename T>
int
ilaenv_LapackTest(int spec, const char *name_, const char *opts,
                  int n1, int n2, int n3, int n4)
{
    using std::complex;

    char name[strlen(name_)+2];
    if (IsSame<T,float>::value) {
        *name = 'S';
    } else if (IsSame<T,double>::value) {
        *name = 'D';
    } else if (IsSame<T,complex<float> >::value) {
        *name = 'C';
    } else if (IsSame<T,complex<double> >::value) {
        *name = 'Z';
    } else {
        ASSERT(0);
    }
    strcpy(name+1, name_);

#if defined CHECK_CXXLAPACK || defined USE_NATIVE_ILAENV
    INTEGER spec_ = spec;
    INTEGER n1_ = n1;
    INTEGER n2_ = n2;
    INTEGER n3_ = n3;
    INTEGER n4_ = n4;
    int result = LAPACK_DECL(ilaenv)(&spec_,
                                     name,
                                     opts,
                                     &n1_,
                                     &n2_,
                                     &n3_,
                                     &n4_,
                                     strlen(name),
                                     strlen(opts));
    return result;
#else
    ASSERT(0);
#endif
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename T>
int
ilaenv(int spec, const char *name, const char *opts,
       int n1, int n2, int n3, int n4)
{
    int info;
#if defined USE_ILAENV_WITH_UNDERSCORE
    info = ilaenv_LapackTest<T>(spec, name, opts, n1, n2, n3, n4);
#else
    info = generic::ilaenv_impl<T>(spec, name, opts, n1, n2, n3, n4);
#endif
    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_ILAENV_TCC
