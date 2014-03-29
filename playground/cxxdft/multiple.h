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

#ifndef PLAYGROUND_CXXDFT_MULTIPLE_H
#define PLAYGROUND_CXXDFT_MULTIPLE_H 1

#define HAVE_CXXDFT_MULTIPLE 1

#include <complex>
#include <playground/cxxdft/direction.h>

namespace cxxdft {

template <typename IndexType, typename VIN, typename VOUT>
    void
    dft_multiple(IndexType n, IndexType m,
                 const VIN *x, IndexType strideX, IndexType distX,
                 VOUT *y, IndexType strideY, IndexType distY,
                 DFTDirection direction);

#ifdef HAVE_FFTW

#ifdef HAVE_FFTW_FLOAT
template <typename IndexType>
    void
    dft_multiple(IndexType n, IndexType m,
                 std::complex<float> *x, IndexType strideX, IndexType distX,
                 std::complex<float> *y, IndexType strideY, IndexType distY,
                 DFTDirection direction);
#endif // HAVE_FFTW_FLOAT

#ifdef HAVE_FFTW_DOUBLE
template <typename IndexType>
    void
    dft_multiple(IndexType n, IndexType m,
                 std::complex<double> *x, IndexType strideX, IndexType distX,
                 std::complex<double> *y, IndexType strideY, IndexType distY,
                 DFTDirection direction);
#endif // HAVE_FFTW_DOUBLE

#ifdef HAVE_FFTW_LONGDOUBLE
template <typename IndexType>
    void
    dft_multiple(IndexType n, IndexType m,
             std::complex<long double> *x, IndexType strideX, IndexType distX,
             std::complex<long double> *y, IndexType strideY, IndexType distY,
             DFTDirection direction);
#endif // HAVE_FFTW_LONGDOUBLE

#ifdef HAVE_FFTW_QUAD
template <typename IndexType>
    void
    dft_multiple(IndexType n, IndexType m,
             std::complex<__float128> *x, IndexType strideX, IndexType distX,
             std::complex<__float128> *y, IndexType strideY, IndexType distY,
             DFTDirection direction);
#endif // HAVE_FFTW_QUAD

#endif

} // namespace cxxfft

#endif // PLAYGROUND_CXXDFT_MULTIPLE_H
