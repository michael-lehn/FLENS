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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_CLASSES_SSE_H
#define PLAYGROUND_CXXBLAS_INTRINSICS_CLASSES_SSE_H 1

#include <playground/cxxblas/intrinsics/includes.h>
#include <playground/cxxblas/intrinsics/classes/functions/functions.h>

#ifdef HAVE_SSE

template <>
class Intrinsics<float, IntrinsicsLevel::SSE> {

public:
	typedef float                            DataType;
	typedef float                            PrimitiveDataType;
	typedef __m128                           IntrinsicsDataType;
	static  const int                        numElements = 4;

	Intrinsics(void)                         {}
	Intrinsics(__m128 val)                   {v = val;}
	Intrinsics(float *a)                     {this->load(a);}
	Intrinsics(float a)                      {this->fill(a);}

	void operator=(float *a)                 {this->load(a);}
	void operator=(float a)                  {this->fill(a);}
    void operator=(__m128 a)                 {v = a; }

	__m128 get(void) const                   {return v;}

    void fill(float a)                       { v = _mm_load1_ps(&a);}
    void load(const float *a)                { v = _mm_load_ps(a);}
    void loadu(const float *a)               { v = _mm_loadu_ps(a);}
	void setZero()                           { v = _mm_setzero_ps();}
    void store(float *a)                     { _mm_store_ps(a, v); }
    void storeu(float *a)                    { _mm_storeu_ps(a, v); }
    void stream(float *a)                    { _mm_stream_ps(a, v); }

private:
	__m128                                   v;

};

template <>
class Intrinsics<double, IntrinsicsLevel::SSE> {

public:
	typedef double                           DataType;
	typedef double                           PrimitiveDataType;
	typedef __m128d                          IntrinsicsDataType;
	static  const int                        numElements = 2;

	Intrinsics(void)                         {}
	Intrinsics(__m128d val)                  {v = val;}
	Intrinsics(double *a)                    {this->load(a);}
	Intrinsics(double a)                     {this->fill(a);}

	void operator=(double *a)                {this->load(a);}
	void operator=(double a)                 {this->fill(a);}
    void operator=(__m128d a)                {v = a; }


	__m128d get(void) const                  {return v;}

    void fill(double a)                      { v = _mm_load1_pd(&a);}
    void load(const double *a)               { v = _mm_load_pd(a);}
    void loadu(const double *a)              { v = _mm_loadu_pd(a);}
	void setZero()                           { v = _mm_setzero_pd();}
    void store(double *a)                    { _mm_store_pd(a, v); }
    void storeu(double *a)                   { _mm_storeu_pd(a, v); }
    void stream(double *a)                   { _mm_stream_pd(a, v); }
private:
	__m128d                                  v;

};

template <>
class Intrinsics<std::complex<float>, IntrinsicsLevel::SSE> {

public:
	typedef std::complex<float>              DataType;
	typedef float                            PrimitiveDataType;
	typedef __m128                           IntrinsicsDataType;
	static  const int                        numElements = 2;

	Intrinsics(void)                         {}
	Intrinsics(__m128 val)                   {v = val;}
	Intrinsics(std::complex<float> *a)       {this->load(a);}
	Intrinsics(float a)                     {this->fill(a);}

	void operator=(std::complex<float> *a)   {this->load(a);}
    void operator=(__m128 a)                 {v = a; }

	__m128 get(void) const                   {return v;}

    void fill(float a)                       { v = _mm_load1_ps(&a); }
    void load(const std::complex<float> *a)  { v = _mm_load_ps(reinterpret_cast<const float* >(a));}
    void loadu(const std::complex<float> *a) { v = _mm_loadu_ps(reinterpret_cast<const float* >(a));}
	void setZero()                           { v = _mm_setzero_ps();}
    void store(std::complex<float> *a)       { _mm_store_ps(reinterpret_cast<float*>(a), v); }
    void storeu(std::complex<float> *a)      { _mm_storeu_ps(reinterpret_cast<float*>(a), v); }
    void stream(std::complex<float> *a)      { _mm_stream_ps(reinterpret_cast<float*>(a), v); }
private:
	__m128                                   v;

};

template <>
class Intrinsics<std::complex<double>, IntrinsicsLevel::SSE> {

public:
	typedef std::complex<double>             DataType;
	typedef double                           PrimitiveDataType;
	typedef __m128d                          IntrinsicsDataType;
	static  const int                        numElements = 1;

	Intrinsics(void)                         {}
	Intrinsics(__m128d val)                  {v = val;}
	Intrinsics(std::complex<double> *a)      {this->load(a);}
    Intrinsics(double a)                     {this->fill(a);}

	void operator=(std::complex<double> *a)  {this->load(a);}
    void operator=(__m128d a)                {v = a; }


	__m128d get(void) const                  {return v;}

    void fill(double a)                      { v = _mm_load1_pd(&a); }
    void load(const std::complex<double> *a) { v = _mm_load_pd(reinterpret_cast<const double* >(a));}
    void loadu(const std::complex<double> *a){ v = _mm_loadu_pd(reinterpret_cast<const double* >(a));}
	void setZero()                           { v = _mm_setzero_pd();}
    void store(std::complex<double> *a)      { _mm_store_pd(reinterpret_cast<double*>(a), v); }
    void storeu(std::complex<double> *a)     { _mm_storeu_pd(reinterpret_cast<double*>(a), v); }
    void stream(std::complex<double> *a)     { _mm_stream_pd(reinterpret_cast<double* >(a), v); }
private:
	__m128d                                  v;

};

#endif // HAVE_SSE

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_CLASSES_CLASSES_H
