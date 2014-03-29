/*
 *   Copyright (c) 2012, Klaus Pototzky
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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_INCLUDES_H
#define PLAYGROUND_CXXBLAS_INTRINSICS_INCLUDES_H 1

enum IntrinsicsLevel {
    AVX  = 2,
    SSE  = 1,
    NONE = 0
};

#ifdef WITH_AVX
#   ifndef HAVE_AVX
#       define HAVE_AVX
#   endif
#   ifndef WITH_SSE
#       define WITH_SSE
#   endif
#   ifndef EMULATE_AVX
#       include <immintrin.h>
#   else
#       include "avxintrin-emu.h"
#   endif
#   ifndef DEFAULT_ALIGNMENT_VALUE
#       define DEFAULT_ALIGNMENT_VALUE     32u
#   endif
#   ifndef DEFAULT_INTRINSIC_LEVEL
#       define DEFAULT_INTRINSIC_LEVEL IntrinsicsLevel::AVX
#   endif
#   ifndef USE_INTRINSIC
#       define USE_INTRINSIC
#   endif
#   ifndef INTRINSIC_NAME
#       define INTRINSIC_NAME        "AVX"
#   endif
#   ifndef _mm256_moveldup_pd
#       define _mm256_moveldup_pd(a) _mm256_permute_pd(a, 0)
#   endif
#   ifndef _mm256_movehdup_pd
#       define _mm256_movehdup_pd(a) _mm256_permute_pd(a, 15)
#   endif
#   ifndef _mm256_load1_ps
#       define _mm256_load1_ps(x) _mm256_set_ps(*x, *x, *x, *x, *x, *x, *x, *x)
#   endif
#   ifndef _mm256_load1_pd
#       define _mm256_load1_pd(x) _mm256_set_pd(*x, *x, *x, *x)
#   endif

#endif

#ifdef WITH_SSE
#   ifndef HAVE_SSE
#       define HAVE_SSE
#   endif
#   ifndef WITH_MMX
#       define WITH_MMX
#   endif
#   include <mmintrin.h>
#   include <xmmintrin.h>
#   include <emmintrin.h>
#   include <pmmintrin.h>
#   include <tmmintrin.h>
#   ifndef DEFAULT_ALIGNMENT_VALUE
#       define DEFAULT_ALIGNMENT_VALUE     16u
#   endif
#   ifndef DEFAULT_INTRINSIC_LEVEL
#       define DEFAULT_INTRINSIC_LEVEL IntrinsicsLevel::SSE
#   endif
#   ifndef USE_INTRINSIC
#       define USE_INTRINSIC
#   endif
#   ifndef INTRINSIC_NAME
#       define INTRINSIC_NAME        "SSE"
#   endif
#   ifndef _mm_moveldup_pd
#       define _mm_moveldup_pd(a) _mm_permute_pd(a, 0)
#   endif
#   ifndef _mm_movehdup_pd
#       define _mm_movehdup_pd(a) _mm_permute_pd(a, 3)
#   endif
#   ifndef _mm_permute_ps
#       define _mm_permute_ps(a,b) _mm_shuffle_ps(a,a,b)
#   endif
#   ifndef _mm_permute_pd
#       define _mm_permute_pd(a,b) _mm_shuffle_pd(a,a,b)
#   endif
#endif

#ifndef DEFAULT_INTRINSIC_LEVEL
#   define DEFAULT_INTRINSIC_LEVEL IntrinsicsLevel::NONE
#endif

#ifndef ASM_COMMENT
#    define ASM_COMMENT(X)  asm("#" X)
#endif

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_INCLUDES_H
