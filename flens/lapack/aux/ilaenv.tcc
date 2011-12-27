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

#ifndef FLENS_LAPACK_AUX_ILAENV_TCC
#define FLENS_LAPACK_AUX_ILAENV_TCC 1

#include <complex>
#include <string>

#include <flens/aux/issame.h>
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

template <typename T>
struct IsNotComplex
{
    static const bool value = true;
};

template <typename T>
struct IsNotComplex<std::complex<T> >
{
    static const bool value = false;
};

template <typename T>
struct IsComplex
{
    static const bool value = !IsNotComplex<T>::value;
};

template <typename T>
int
ilaenv_generic(int spec, const char *name, const char *opts,
               int n1, int n2, int n3, int n4)
{
    int result = -1;
    const char *c2 = name;
    const char *c3 = name + 2;
    const char *c4 = name + 1;

    int nb, nbMin, nx;

    switch (spec) {
//
//      optimal blocksize
//
        case 1:
            nb = -1;
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
            nx = 0;
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
            result = nx;
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

//== interface for native lapack ===============================================

#if defined CHECK_CXXLAPACK || defined USE_NATIVE_ILAENV

#ifndef LAPACK_DECL
#   define  LAPACK_DECL(x)    x##_

#   define  INTEGER           int

extern "C" {

INTEGER
LAPACK_DECL(ilaenv)(const INTEGER   *SPEC,
                    const char      *NAME,
                    const char      *OPTS,
                    const INTEGER   *N1,
                    const INTEGER   *N2,
                    const INTEGER   *N3,
                    const INTEGER   *N4,
                    int             NAME_LEN,
                    int             OPTS_LEN);

} // extern "C"

#endif


template <typename T>
int
ilaenv_LapackTest(int spec, const char *_name, const char *_opts,
                  int n1, int n2, int n3, int n4)
{
    using std::string;
    using std::complex;

    string opts(_opts);
    string name;
    if (IsSame<T,float>::value) {
        name = string("S") + string(_name);
    } else if (IsSame<T,double>::value) {
        name = string("D") + string(_name);
    } else if (IsSame<T,complex<float> >::value) {
        name = string("C") + string(_name);
    } else if (IsSame<T,complex<double> >::value) {
        name = string("Z") + string(_name);
    } else {
        ASSERT(0);
    }

#if defined CHECK_CXXLAPACK || defined USE_NATIVE_ILAENV
    int result = LAPACK_DECL(ilaenv)(&spec, name.c_str(), opts.c_str(),
                                     &n1, &n2, &n3, &n4,
                                     strlen(name.c_str()),
                                     strlen(opts.c_str()));

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
#if defined CHECK_CXXLAPACK || defined USE_NATIVE_ILAENV
    return ilaenv_LapackTest<T>(spec, name, opts, n1, n2, n3, n4);
#else
    return ilaenv_generic<T>(spec, name, opts, n1, n2, n3, n4);
#endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_ILAENV_TCC
