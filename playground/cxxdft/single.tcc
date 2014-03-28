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
 */

#ifndef PLAYGROUND_CXXDFT_SINGLE_TCC
#define PLAYGROUND_CXXDFT_SINGLE_TCC 1

#include <cmath>
#include <string.h>
#include <cxxblas/cxxblas.h>
#include <flens/auxiliary/auxiliary.h>
#include <playground/cxxdft/direction.h>

namespace cxxdft {

template <typename T>
typename flens::RestrictTo<flens::IsInteger<T>::value,
                           bool>::Type
isPowerOfTwo (T x)
{
  return ((x != 0) && ((x & (~x + 1)) == x));
};

template<typename IndexType, typename VIN, typename VOUT>
void
fft_single_generic(IndexType N,
                   const VIN *x, IndexType incX,
                   VOUT *y, IndexType incY,
                   DFTDirection direction)
{
    CXXBLAS_DEBUG_OUT("fft_single_generic");

    typedef typename flens::ComplexTrait<VOUT>::PrimitiveType PT;

    if ( N <= 1 ) {
        return;
    }

    const PT        factor = (direction==DFTDirection::Forward ? PT(-2*M_PI/N) : PT(2*M_PI/N) );
    const PT        one(1);
    const IndexType Nhalf = N/2;
    if (N == 1) {
        y[0] = x[0];
        return;
    }

    VOUT *even = new VOUT[Nhalf];
    VOUT *odd  = new VOUT[Nhalf];

    for(IndexType k=0, iX=0; k < Nhalf; ++k, iX+=2*incX) {
        even[k] = x[iX];
        odd[k]  = x[iX+incX];
    }

    fft_single_generic(Nhalf, even, 1, even, 1, direction);
    fft_single_generic(Nhalf, odd , 1, odd , 1, direction);

    for (IndexType k = 0, iY=0; k < Nhalf; ++k, iY+=incY)
    {
        VOUT t = std::polar(one, factor * k) * odd[k];
        y[iY           ] = even[k] + t;
        y[iY+incY*Nhalf] = even[k] - t;
    }

    delete[] even;
    delete[] odd;

    return;
}

template<typename IndexType, typename VIN, typename VOUT>
void
dft_single_generic(IndexType N,
                   const VIN *x, IndexType incX,
                   VOUT *y, IndexType incY,
                   DFTDirection direction)
{

    CXXBLAS_DEBUG_OUT("dft_single_generic");

    //
    // Use Cooley-Tukey algorithm if possible
    //
    if ( isPowerOfTwo(N) ) {

        fft_single_generic(N, x, incX, y, incY, direction);

    } else {

        typedef typename flens::ComplexTrait<VOUT>::PrimitiveType PT;

        const PT factor = (direction==DFTDirection::Forward ? PT(-2*M_PI/N) : PT(2*M_PI/N));

        for (IndexType i=0, iY=0; i<N; ++i, iY+=incY) {
            VOUT tmp(0);
            for (IndexType j=0, iX=0; j<N; ++j, iX+=incX) {
                tmp += std::polar(PT(1),factor*i*j)*x[iX];
            }
            y[iY] = tmp;
        }

    }
}

template <typename IndexType, typename VIN, typename VOUT>
void
dft_single(IndexType n,
           const VIN *x, IndexType incX,
           VOUT *y, IndexType incY,
           DFTDirection direction)
{
    if (incX<0) {
        x -= incX*(n-1);
    }
    if (incY<0) {
        y -= incY*(n-1);
    }
    dft_single_generic(n, x, incX, y, incY, direction);
}

#ifdef HAVE_FFTW

#ifdef HAVE_FFTW_FLOAT
template <typename IndexType>
void
dft_single(IndexType n,
           std::complex<float> *x, IndexType incX,
           std::complex<float> *y, IndexType incY,
           DFTDirection direction)
{
    CXXBLAS_DEBUG_OUT("dft_single [FFTW interface, float]");

    fftwf_plan p = NULL;

#   if defined(FFTW_WISDOM_IMPORT) && !defined(WITH_MKLBLAS)
        ASSERT(strcmp(FFTW_WISDOM_FILENAME,""));
        fftwf_import_wisdom_from_filename(FFTW_WISDOM_FILENAME);
#   endif

    if ( incX==1 && incX==1 ) {
        p = fftwf_plan_dft_1d(n,
                              reinterpret_cast<fftwf_complex*>(x),
                              reinterpret_cast<fftwf_complex*>(y),
                              direction, FFTW_PLANNER_FLAG);
    } else {
        p = fftwf_plan_many_dft(1, &n, 1,
                                reinterpret_cast<fftwf_complex*>(x), NULL, incX, 0,
                                reinterpret_cast<fftwf_complex*>(y), NULL, incY, 0,
                                direction, FFTW_PLANNER_FLAG);
    }

    fftwf_execute(p);

#   if defined(FFTW_WISDOM_EXPORT) && !defined(WITH_MKLBLAS)
        ASSERT(strcmp(FFTW_WISDOM_FILENAME,""));
        fftwf_export_wisdom_to_filename(FFTW_WISDOM_FILENAME);
#   endif

    fftwf_destroy_plan(p);
}

#endif // HAVE_FFTW_FLOAT

#ifdef HAVE_FFTW_DOUBLE
template <typename IndexType>
void
dft_single(IndexType n,
           std::complex<double> *x, IndexType incX,
           std::complex<double> *y, IndexType incY,
           DFTDirection direction)
{
    CXXBLAS_DEBUG_OUT("dft_single [FFTW interface, double]");

    fftw_plan p = NULL;

#   if defined(FFTW_WISDOM_IMPORT) && !defined(WITH_MKLBLAS)
        ASSERT(strcmp(FFTW_WISDOM_FILENAME,""));
        fftw_import_wisdom_from_filename(FFTW_WISDOM_FILENAME);
#   endif

    if ( incX==1 && incX==1 ) {
        p = fftw_plan_dft_1d(n,
                             reinterpret_cast<fftw_complex*>(x),
                             reinterpret_cast<fftw_complex*>(y),
                             direction, FFTW_PLANNER_FLAG);
    } else {
        p = fftw_plan_many_dft(1, &n, 1,
                               reinterpret_cast<fftw_complex*>(x), NULL, incX, 0,
                               reinterpret_cast<fftw_complex*>(y), NULL, incY, 0,
                               direction, FFTW_PLANNER_FLAG);
    }

    fftw_execute(p);

#   if defined(FFTW_WISDOM_EXPORT) && !defined(WITH_MKLBLAS)
        ASSERT(strcmp(FFTW_WISDOM_FILENAME,""));
        fftw_export_wisdom_to_filename(FFTW_WISDOM_FILENAME);
#   endif

    fftw_destroy_plan(p);
}

#endif // HAVE_FFTW_DOUBLE

#ifdef HAVE_FFTW_LONGDOUBLE
template <typename IndexType>
void
dft_single(IndexType n,
           std::complex<long double> *x, IndexType incX,
           std::complex<long double> *y, IndexType incY,
           DFTDirection direction)
{
    CXXBLAS_DEBUG_OUT("dft_single [FFTW interface, long double]");

    fftwl_plan p = NULL;

#   ifdef FFTW_WISDOM_IMPORT
        ASSERT(strcmp(FFTW_WISDOM_FILENAME,""));
        fftwl_import_wisdom_from_filename(FFTW_WISDOM_FILENAME);
#   endif

    if ( incX==1 && incX==1 ) {
        p = fftwl_plan_dft_1d(n,
                              reinterpret_cast<fftwl_complex*>(x),
                              reinterpret_cast<fftwl_complex*>(y),
                              direction, FFTW_PLANNER_FLAG);
    } else {
        p = fftwl_plan_many_dft(1, &n, 1,
                                reinterpret_cast<fftwl_complex*>(x), NULL, incX, 0,
                                reinterpret_cast<fftwl_complex*>(y), NULL, incY, 0,
                                direction, FFTW_PLANNER_FLAG);
    }

    fftwl_execute(p);

#   ifdef FFTW_WISDOM_EXPORT
        ASSERT(strcmp(FFTW_WISDOM_FILENAME,""));
        fftwl_export_wisdom_to_filename(FFTW_WISDOM_FILENAME);
#   endif

    fftwl_destroy_plan(p);

}

#endif // HAVE_FFTW_LONGDOUBLE

#ifdef HAVE_FFTW_QUAD
template <typename IndexType>
void
dft_single(IndexType n,
           std::complex<__float128> *x, IndexType incX,
           std::complex<__float128> *y, IndexType incY,
           DFTDirection direction)
{
    CXXBLAS_DEBUG_OUT("dft_single [FFTW interface, quad]");

    fftwq_plan p = NULL;

#   ifdef FFTW_WISDOM_IMPORT
        ASSERT(strcmp(FFTW_WISDOM_FILENAME,""));
        fftwq_import_wisdom_from_filename(FFTW_WISDOM_FILENAME);
#   endif

    if ( incX==1 && incX==1 ) {
        p = fftwq_plan_dft_1d(n,
                              reinterpret_cast<fftwq_complex*>(x),
                              reinterpret_cast<fftwq_complex*>(y),
                              direction, FFTW_PLANNER_FLAG);
    } else {
        p = fftwq_plan_many_dft(1, &n, 1,
                                reinterpret_cast<fftwq_complex*>(x), NULL, incX, 0,
                                reinterpret_cast<fftwq_complex*>(y), NULL, incY, 0,
                                direction, FFTW_PLANNER_FLAG);
    }

    fftwq_execute(p);

#   ifdef FFTW_WISDOM_EXPORT
        ASSERT(strcmp(FFTW_WISDOM_FILENAME,""));
        fftwq_export_wisdom_to_filename(FFTW_WISDOM_FILENAME);
#   endif

    fftwq_destroy_plan(p);

}

#endif // HAVE_FFTW_QUAD

#endif


} // namespace cxxdft

#endif // PLAYGROUND_CXXDFT_SINGLE_TCC
