/*
 *   Copyright (c) 2013, Klaus Pototzky
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
 *
 *   Source:
 *   http://nickstrupat.blogspot.ca/2012/12/i-wrote-this-function-for-cache-line.html
 *
 */

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_AUXILIARY_CACHESIZE_H
#define PLAYGROUND_CXXBLAS_INTRINSICS_AUXILIARY_CACHESIZE_H 1

#ifdef USE_INTRINSIC

#if defined(__linux__)
#include <unistd.h>

inline
size_t get_l1_cache_size()
{
    static size_t l1_cache = 0;
    if (l1_cache==0) {
        l1_cache=sysconf(SC_LEVEL1_DCACHE_SIZE_);
    }
    return l1_cache;
}

inline
size_t get_l2_cache_size()
{
    static size_t l2_cache = 0;
    if (l2_cache==0) {
        l2_cache=sysconf(SC_LEVEL2_CACHE_SIZE_);
    }
    return l2_cache;
}

inline
size_t get_l3_cache_size()
{
    static size_t l3_cache = 0;
    if (l3_cache==0) {
        l3_cache=sysconf(SC_LEVEL3_CACHE_SIZE_);
    }
    return l3_cache;
}

#elif defined (__APPLE__)

#include <sys/sysctl.h>

inline
size_t get_l1_cache_size()
{
    size_t line_size = 0;
    size_t sizeof_line_size = sizeof(line_size);
    sysctlbyname("hw.l1dcachesize", &line_size, &sizeof_line_size, 0, 0);
    return line_size;
}

inline
size_t get_l2_cache_size()
{
    size_t line_size = 0;
    size_t sizeof_line_size = sizeof(line_size);
    sysctlbyname("hw.l2cachesize", &line_size, &sizeof_line_size, 0, 0);
    return line_size;
}

inline
size_t get_l3_cache_size()
{
    size_t line_size = 0;
    size_t sizeof_line_size = sizeof(line_size);
    sysctlbyname("hw.l3cachesize", &line_size, &sizeof_line_size, 0, 0);
    return line_size;
}

#else
// Cache sizes have to be set manually

inline
size_t get_l1_cache_size()
{
    return L1_CACHE_SIZE;
}

inline
size_t get_l2_cache_size()
{
    return L2_CACHE_SIZE;
}

inline
size_t get_l3_cache_size()
{
    return L3_CACHE_SIZE;
}
#endif

#endif // USE_INTRINSIC

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_AUXILIARY_CACHESIZE_H
