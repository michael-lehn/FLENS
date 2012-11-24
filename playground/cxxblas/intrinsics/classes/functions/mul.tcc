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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_CLASSES_FUNCTIONS_MUL_TCC
#define PLAYGROUND_CXXBLAS_INTRINSICS_CLASSES_FUNCTIONS_MUL_TCC 1

#include <playground/cxxblas/intrinsics/includes.h>

#ifdef HAVE_SSE

//--- Mul
Intrinsics<float, IntrinsicsLevel::SSE>
inline _intrinsic_mul(const Intrinsics<float, IntrinsicsLevel::SSE> &x,
                      const Intrinsics<float, IntrinsicsLevel::SSE> &y)
{
	return Intrinsics<float, IntrinsicsLevel::SSE>(_mm_mul_ps(x.get(), y.get()));
}

Intrinsics<double, IntrinsicsLevel::SSE>
inline _intrinsic_mul(const Intrinsics<double, IntrinsicsLevel::SSE> &x,
                      const Intrinsics<double, IntrinsicsLevel::SSE> &y)
{
	return Intrinsics<double, IntrinsicsLevel::SSE>(_mm_mul_pd(x.get(), y.get()));
}

Intrinsics<std::complex<float>, IntrinsicsLevel::SSE>
inline _intrinsic_mul(const Intrinsics<std::complex<float>, IntrinsicsLevel::SSE> &x,
                      const Intrinsics<float, IntrinsicsLevel::SSE> &y)
{
	return Intrinsics<std::complex<float>, IntrinsicsLevel::SSE>(_mm_mul_ps(x.get(), y.get()));
}

Intrinsics<std::complex<float>, IntrinsicsLevel::SSE>
inline _intrinsic_mul(const Intrinsics<float, IntrinsicsLevel::SSE> &x,
                      const Intrinsics<std::complex<float>, IntrinsicsLevel::SSE> &y)
{
	return Intrinsics<std::complex<float>, IntrinsicsLevel::SSE>(_mm_mul_ps(x.get(), y.get()));
}

Intrinsics<std::complex<double>, IntrinsicsLevel::SSE>
inline _intrinsic_mul(const Intrinsics<std::complex<double>, IntrinsicsLevel::SSE> &x,
                      const Intrinsics<double, IntrinsicsLevel::SSE> &y)
{
	return Intrinsics<std::complex<double>, IntrinsicsLevel::SSE>(_mm_mul_pd(x.get(), y.get()));
}

Intrinsics<std::complex<double>, IntrinsicsLevel::SSE>
inline _intrinsic_mul(const Intrinsics<double, IntrinsicsLevel::SSE> &x,
                      const Intrinsics<std::complex<double>, IntrinsicsLevel::SSE> &y)
{
	return Intrinsics<std::complex<double>, IntrinsicsLevel::SSE>(_mm_mul_pd(x.get(), y.get()));
}

#endif // HAVE_SSE


#ifdef HAVE_AVX

//--- Mul
Intrinsics<float, IntrinsicsLevel::AVX>
inline _intrinsic_mul(const Intrinsics<float, IntrinsicsLevel::AVX> &x,
		              const Intrinsics<float, IntrinsicsLevel::AVX> &y)
{
	return Intrinsics<float, IntrinsicsLevel::AVX>(_mm256_mul_ps(x.get(), y.get()));
}

Intrinsics<double, IntrinsicsLevel::AVX>
inline _intrinsic_mul(const Intrinsics<double, IntrinsicsLevel::AVX> &x,
		              const Intrinsics<double, IntrinsicsLevel::AVX> &y)
{
	return Intrinsics<double, IntrinsicsLevel::AVX>(_mm256_mul_pd(x.get(), y.get()));
}

Intrinsics<std::complex<float>, IntrinsicsLevel::AVX>
inline _intrinsic_mul(const Intrinsics<std::complex<float>, IntrinsicsLevel::AVX> &x,
		              const Intrinsics<float, IntrinsicsLevel::AVX> &y)
{
	return Intrinsics<std::complex<float>, IntrinsicsLevel::AVX>(_mm256_mul_ps(x.get(), y.get()));
}

Intrinsics<std::complex<float>, IntrinsicsLevel::AVX>
inline _intrinsic_mul(const Intrinsics<float, IntrinsicsLevel::AVX> &x,
		              const Intrinsics<std::complex<float>, IntrinsicsLevel::AVX> &y)
{
	return Intrinsics<std::complex<float>, IntrinsicsLevel::AVX>(_mm256_mul_ps(x.get(), y.get()));
}

Intrinsics<std::complex<double>, IntrinsicsLevel::AVX>
inline _intrinsic_mul(const Intrinsics<std::complex<double>, IntrinsicsLevel::AVX> &x,
		              const Intrinsics<double, IntrinsicsLevel::AVX> &y)
{
	return Intrinsics<std::complex<double>, IntrinsicsLevel::AVX>(_mm256_mul_pd(x.get(), y.get()));
}

Intrinsics<std::complex<double>, IntrinsicsLevel::AVX>
inline _intrinsic_mul(const Intrinsics<double, IntrinsicsLevel::AVX> &x,
		             const Intrinsics<std::complex<double>, IntrinsicsLevel::AVX> &y)
{
	return Intrinsics<std::complex<double>, IntrinsicsLevel::AVX>(_mm256_mul_pd(x.get(), y.get()));
}

#endif // HAVE_AVX

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_CLASSES_FUNCTIONS_MUL_TCC
