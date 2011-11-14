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

template <typename T>
int
ilaenv_generic(int spec, const char *_name, const char *_opts,
               int n1, int n2, int n3, int n4)
{
    std::cerr << "spec = " << spec
              << ", name = " << _name
              << ", opts = " << _opts
              << ", n1 = " << n1
              << ", n2 = " << n2
              << ", n3 = " << n3
              << ", n4 = " << n4
              << std::endl;
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

#   ifdef LAPACK_DECL
    int result = LAPACK_DECL(ilaenv)(&spec, name.c_str(), opts.c_str(),
                                     &n1, &n2, &n3, &n4,
                                     strlen(name.c_str()),
                                     strlen(opts.c_str()));
    std::cerr << "ilaenv: result = " << result << std::endl;
#   else
    int result = 1;
#   endif

    return result;
}

//== interface for native lapack ===============================================

#if defined CHECK_CXXLAPACK || defined USE_NATIVE_ILAENV

template <typename T>
int
ilaenv_native(int spec, const char *_name, const char *_opts,
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
#else
    ASSERT(0);
#endif

    //std::cerr << "ilaenv: result = " << result << std::endl;

    return result;
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename T>
int
ilaenv(int spec, const char *name, const char *opts,
       int n1, int n2, int n3, int n4)
{
#if defined CHECK_CXXLAPACK || defined USE_NATIVE_ILAENV
    return ilaenv_native<T>(spec, name, opts, n1, n2, n3, n4);
#else
    return ilaenv_generic<T>(spec, name, opts, n1, n2, n3, n4);
#endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_ILAENV_TCC
