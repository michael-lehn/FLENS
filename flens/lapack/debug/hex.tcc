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

#ifndef FLENS_LAPACK_DEBUG_HEX_TCC
#define FLENS_LAPACK_DEBUG_HEX_TCC 1

#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <flens/auxiliary/auxiliary.h>

namespace flens { namespace lapack {

template <typename T>
std::string
hex(const T &x)
{
    std::stringstream out;

    const unsigned int n = sizeof(T);
    const unsigned char *c = reinterpret_cast<const unsigned char*>(&x);

    for (unsigned int i=0; i<n; ++i) {
//
//      Lehn: No idea why I have to copy c[i] into 'tmp' ...
//
        unsigned int tmp = c[i];
        out << "0x";
        out.width(2);
        out.fill('0');
        out << std::hex << tmp;
        if (i<n-1) {
            out << " ";
        }
    }
    return out.str();
}

template <typename T>
T
hex(const char *str)
{
    const unsigned int n = sizeof(T);
    ASSERT(strlen(str)==n*5-1);

    T result;
    unsigned char *c = reinterpret_cast<unsigned char*>(&result);

    for (unsigned int i=0; i<n; ++i, str+=5) {
        c[i] = strtol(str, 0, 16);
    }

    return result;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_DEBUG_HEX_TCC
