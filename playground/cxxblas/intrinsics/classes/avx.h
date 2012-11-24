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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_CLASSES_AVX_H
#define PLAYGROUND_CXXBLAS_INTRINSICS_CLASSES_AVX_H 1

#include <playground/cxxblas/intrinsics/includes.h>
#include <playground/cxxblas/intrinsics/classes/functions/functions.h>

#ifdef HAVE_AVX

template <>
class Intrinsics<float, IntrinsicsLevel::AVX> {

public:
	typedef float                            DataType;
	typedef float                            PrimitiveDataType;
	typedef __m256                           IntrinsicsDataType;
	static  const int                        numElements = 8;

	Intrinsics(void)                         {}
	Intrinsics(__m256 val)                   {v = val;}
	Intrinsics(float *a)                     {this->load(a);}
	Intrinsics(float a)                      {this->fill(a);}

	void operator=(float *a)                 {this->load(a);}
	void operator=(float a)                  {this->fill(a);}
    void operator=(__m256 a)                 {v = a; }

	__m256 get(void) const                   {return v;}

    void fill(float a)                       { v = _mm256_broadcast_ss(&a); }
    void load(const float *a)                { v = _mm256_load_ps(a);  }
    void loadu(const float *a)               { v = _mm256_loadu_ps(a);  }
	void setZero()                           { v = _mm256_setzero_ps(); }
    void store(float *a)                     { _mm256_store_ps(a, v);  }
    void storeu(float *a)                    { _mm256_storeu_ps(a, v);  }
    void stream(float *a)                    { _mm256_stream_ps(a, v); }

private:
	__m256                                   v;

};

template <>
class Intrinsics<double, IntrinsicsLevel::AVX> {

public:
	typedef double                           DataType;
	typedef double                           PrimitiveDataType;
	typedef __m256d                          IntrinsicsDataType;
	static  const int                        numElements = 4;

	Intrinsics(void)                         {}
	Intrinsics(__m256d val)                  {v = val;}
	Intrinsics(double *a)                    {this->load(a);}
	Intrinsics(double a)                     {this->fill(a);}

	void operator=(double *a)                {this->load(a);}
	void operator=(double a)                 {this->fill(a);}
    void operator=(__m256d a)                {v = a; }


	__m256d get(void) const                  {return v;}

    void fill(double a)                      { v = _mm256_broadcast_sd(&a); }
    void load(const double *a)               { v = _mm256_load_pd(a);  }
    void loadu(const double *a)              { v = _mm256_loadu_pd(a);  }
	void setZero()                           { v = _mm256_setzero_pd(); }
    void store(double *a)                    { _mm256_store_pd(a, v);  }
    void storeu(double *a)                   { _mm256_storeu_pd(a, v);  }
    void stream(double *a)                   { _mm256_stream_pd(a, v); }

private:
	__m256d                                  v;

};

template <>
class Intrinsics<std::complex<float>, IntrinsicsLevel::AVX> {

public:
	typedef std::complex<float>              DataType;
	typedef float                            PrimitiveDataType;
	typedef __m256                           IntrinsicsDataType;
	static  const int                        numElements = 4;

	Intrinsics(void)                         {}
	Intrinsics(__m256 val)                   {v = val;}
	Intrinsics(std::complex<float> *a)       {this->load(a);}
    Intrinsics(float a)                      {this->fill(a);}

	void operator=(std::complex<float> *a)   {this->load(a);}
    void operator=(__m256 a)                 {v = a; }

	__m256 get(void) const                   {return v;}

    void fill(float a)                       { v = _mm256_broadcast_ss(&a); }
    void load(const std::complex<float> *a)  { v = _mm256_load_ps(reinterpret_cast<const float* >(a));}
    void loadu(const std::complex<float> *a) { v = _mm256_loadu_ps(reinterpret_cast<const float* >(a));}
	void setZero()                           { v = _mm256_setzero_ps();}
    void store(std::complex<float> *a)       { _mm256_store_ps(reinterpret_cast<float* >(a), v); }
    void storeu(std::complex<float> *a)      { _mm256_storeu_ps(reinterpret_cast<float* >(a), v); }
    void stream(std::complex<float> *a)      { _mm256_stream_ps(reinterpret_cast<float *>(a), v); }

private:
	__m256                                   v;

};

template <>
class Intrinsics<std::complex<double>, IntrinsicsLevel::AVX> {

public:
	typedef std::complex<double>             DataType;
	typedef double                           PrimitiveDataType;
	typedef __m256d                          IntrinsicsDataType;
	static  const int                        numElements = 2;

	Intrinsics(void)                         {}
	Intrinsics(__m256d val)                  {v = val;}
	Intrinsics(std::complex<double> *a)      {this->load(a);}
    Intrinsics(double a)                     {this->fill(a);}

	void operator=(std::complex<double> *a)  {this->load(a);}
    void operator=(__m256d a)                {v = a; }


	__m256d get(void) const                  {return v;}

    void fill(double a)                      { v = _mm256_broadcast_sd(&a); }
    void load(const std::complex<double> *a) { v = _mm256_load_pd(reinterpret_cast<const double* >(a));}
    void loadu(const std::complex<double> *a){ v = _mm256_loadu_pd(reinterpret_cast<const double* >(a));}
	void setZero()                           { v = _mm256_setzero_pd();}
    void store(std::complex<double> *a)      { _mm256_storeu_pd(reinterpret_cast<double* >(a), v);}
    void storeu(std::complex<double> *a)     { _mm256_storeu_pd(reinterpret_cast<double* >(a), v);}
    void stream(std::complex<double> *a)     { _mm256_stream_pd(reinterpret_cast<double* >(a), v); }
private:
	__m256d                                  v;

};

#endif // HAVE_AVX

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_CLASSES_AVX_H
