/*
 *   Copyright (c) 2011, Iris Haecker
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


/*
 *   Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.
 *   
 *   $COPYRIGHT$
 *   
 *   Additional copyrights may follow
 *   
 *   $HEADER$
 *   
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions are
 *   met:
 *   
 *   - Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer. 
 *     
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer listed
 *     in this license in the documentation and/or other materials
 *     provided with the distribution.
 *     
 *   - Neither the name of the copyright holders nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
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

#ifndef CXXLAPACK_AUX_LARTG_TCC
#define CXXLAPACK_AUX_LARTG_TCC 1

#include <cxxblas/cxxblas.h>
#include <cxxlapack/cxxlapack.h>

namespace cxxlapack {

template <typename T>
void
lartg(const T &f, const T &g, T &c, T &s, T &r)
{
    typedef typename LartgInt<T>::Type Integer;

    static const T safeMin = lamch<T>(SafeMin);
    static const T eps = lamch<T>(Eps);
    static const T safeMin2 = pow(lamch<T>(Base),
                                  Integer(log(safeMin/eps)
                                            / log(lamch<T>(Base))
                                            / T(2))
                                 );
    static const T safeMax2 = T(1) / safeMin2;

    if (g==T(0)) {
        c = T(1);
        s = T(0);
        r = f;
    } else if (f==T(0)) {
        c = T(0);
        s = T(1);
        r = g;
    } else {
        T f_ = f;
        T g_ = g;
        T scale = max(abs(f_), abs(g_));
        if (scale>=safeMax2) {
            int count = 0;
            do {
                ++count;
                f_ *= safeMin2;
                g_ *= safeMin2;
                scale = max(abs(f_), abs(g_));
            } while (scale>=safeMax2);
            r = sqrt(f_*f_ + g_*g_);
            c = f_ / r;
            s = g_ / r;
            for (int i=0; i<count; ++i) {
                r *= safeMax2;
            }
        } else if (scale<=safeMin2) {
            int count = 0;
            do {
                ++count;
                f_ *= safeMax2;
                g_ *= safeMax2;
                scale = max(abs(f_), abs(g_));
            } while (scale<=safeMin2);
            r = sqrt(f_*f_ + g_*g_);
            c = f_ / r;
            s = g_ / r;
            for (int i=0; i<count; ++i) {
                r *= safeMin2;
            }
        } else {
            r = sqrt(f_*f_ + g_*g_);
            c = f_ / r;
            s = g_ / r;
        }
        if ((abs(f)>abs(g)) && (c<T(0))) {
            c = -c;
            s = -s;
            r = -r;
        }
    }
}

} // namespace cxxlapack

#endif // CXXLAPACK_AUX_LARTG_TCC
