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
       INTEGER FUNCTION IPARMQ( spec, NAME, OPTS, N, ILO, IHI, LWORK )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_IPARMQ_TCC
#define FLENS_LAPACK_IMPL_IPARMQ_TCC 1

#include <complex>
#include <string>

#include <flens/auxiliary/issame.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== publich interface / generic lapack implementation =========================

template <typename T>
int
iparmq(int spec, const char *, const char *, int, int iLo, int iHi, int)
{
    using std::log;
    using std::max;

    const int INMIN = 12;
    const int INWIN = 13;
    const int INIBL = 14;
    const int ISHFTS = 15;
    const int IACC22 = 16;

    const float Two(2);

    const int nMin = 75;
    const int k22Min = 14;
    const int kacMin = 14;
    const int nibble = 14;
    const int knwSwp = 500;

    int nh = -1,
        ns = -1,
        result = -1;

    if ((spec==ISHFTS) || (spec==INWIN) || (spec==IACC22)) {
//
//      ==== Set the number simultaneous shifts ====
//
        nh = iHi - iLo + 1;
        ns = 2;
        if (nh>=30) {
            ns = 4;
        }
        if (nh>=60) {
            ns = 10;
        }
        if (nh>=150) {
            ns = max(10, int(nh / nint(log(float(nh)) / log(Two))));
        }
        if (nh>=590) {
            ns = 64;
        }
        if (nh>=3000) {
            ns = 128;
        }
        if (nh>=6000) {
            ns = 256;
        }
        ns = max(2, ns - (ns % 2));
    }

    if (spec==INMIN) {
//
//
//      ===== Matrices of order smaller than NMIN get sent
//      .     to xLAHQR, the classic double shift algorithm.
//      .     This must be at least 11. ====
//
        result = nMin;

    } else if (spec==INIBL) {
//
//      ==== INIBL: skip a multi-shift qr iteration and
//      .    whenever aggressive early deflation finds
//      .    at least (NIBBLE*(window size)/100) deflations. ====
//
        result = nibble;

    } else if (spec==ISHFTS) {
//
//      ==== NSHFTS: The number of simultaneous shifts =====
//
        result = ns;

    } else if (spec==INWIN) {
//
//      ==== NW: deflation window size.  ====
//
        if (nh<=knwSwp) {
            result = ns;
        } else {
            result = 3*ns / 2;
        }

    } else if (spec==IACC22) {
//
//      ==== IACC22: Whether to accumulate reflections
//      .     before updating the far-from-diagonal elements
//      .     and whether to use 2-by-2 block structure while
//      .     doing it.  A small amount of work could be saved
//      .     by making this choice dependent also upon the
//      .     NH=IHI-ILO+1.
//
        result = 0;
        if (ns>=kacMin) {
            result = 1;
        }
        if (ns>=k22Min) {
            result = 2;
        }

    } else {
//
//      ===== invalid value of ispec =====
//
        result = -1;

    }

    return result;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_IPARMQ_TCC
