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

#ifndef PLAYGROUND_CXXDFT_MULTIPLE_TCC
#define PLAYGROUND_CXXDFT_MULTIPLE_TCC 1

#include <cmath>
#include <string.h>
#include <flens/auxiliary/auxiliary.h>
#include <playground/cxxdft/single.tcc>
#include <playground/cxxdft/direction.h>

namespace cxxdft {

template <typename IndexType, typename VIN, typename VOUT>
void
dft_multiple(IndexType n, IndexType m,
             const VIN *x, IndexType strideX, IndexType distX,
             VOUT *y, IndexType strideY, IndexType distY,
             DFTDirection direction)
{
    CXXBLAS_DEBUG_OUT("dft_multiple");

    for (IndexType i=0; i<m; ++i) {
        dft_single_generic(n, x+i*distX, strideX, y+i*distY, strideY,
                           direction);
    }

}

#ifdef HAVE_FFTW

#ifdef HAVE_FFTW_FLOAT

template <typename IndexType>
void
dft_multiple(IndexType n, IndexType m,
             std::complex<float> *x, IndexType strideX, IndexType distX,
             std::complex<float> *y, IndexType strideY, IndexType distY,
             DFTDirection direction)
{
    CXXBLAS_DEBUG_OUT("dft_multiple [FFTW interface, double]");

#   if defined(FFTW_WISDOM_IMPORT) && !defined(WITH_MKLBLAS)
        ASSERT(strcmp(FFTW_WISDOM_FILENAME,""));
        fftwf_import_wisdom_from_filename(FFTW_WISDOM_FILENAME);
#   endif

    fftwf_plan p = fftwf_plan_many_dft(1, &n, m,
                                      reinterpret_cast<fftwf_complex*>(x), NULL, strideX, distX,
                                      reinterpret_cast<fftwf_complex*>(y), NULL, strideY, distY,
                                      direction, FFTW_PLANNER_FLAG);
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
dft_multiple(IndexType n, IndexType m,
             std::complex<double> *x, IndexType strideX, IndexType distX,
             std::complex<double> *y, IndexType strideY, IndexType distY,
             DFTDirection direction)
{
    CXXBLAS_DEBUG_OUT("dft_multiple [FFTW interface, double]");

#   if defined(FFTW_WISDOM_IMPORT) && !defined(WITH_MKLBLAS)
        ASSERT(strcmp(FFTW_WISDOM_FILENAME,""));
        fftw_import_wisdom_from_filename(FFTW_WISDOM_FILENAME);
#   endif

    fftw_plan p = fftw_plan_many_dft(1, &n, m,
                                     reinterpret_cast<fftw_complex*>(x), NULL, strideX, distX,
                                     reinterpret_cast<fftw_complex*>(y), NULL, strideY, distY,
                                     direction, FFTW_PLANNER_FLAG);
    fftw_execute(p);

#  if defined(FFTW_WISDOM_EXPORT) && !defined(WITH_MKLBLAS)
        ASSERT(strcmp(FFTW_WISDOM_FILENAME,""));
        fftw_export_wisdom_to_filename(FFTW_WISDOM_FILENAME);
#   endif

    fftw_destroy_plan(p);
}

#endif // HAVE_FFTW_DOUBLE

#ifdef HAVE_FFTW_LONGDOUBLE

template <typename IndexType>
void
dft_multiple(IndexType n, IndexType m,
             std::complex<long double> *x, IndexType strideX, IndexType distX,
             std::complex<long double> *y, IndexType strideY, IndexType distY,
             DFTDirection direction)
{
    CXXBLAS_DEBUG_OUT("dft_multiple [FFTW interface, long double]");

#   ifdef FFTW_WISDOM_IMPORT
        ASSERT(strcmp(FFTW_WISDOM_FILENAME,""));
        fftwl_import_wisdom_from_filename(FFTW_WISDOM_FILENAME);
#   endif

    fftwl_plan p = fftwl_plan_many_dft(1, &n, m,
                                      reinterpret_cast<fftwl_complex*>(x), NULL, strideX, distX,
                                      reinterpret_cast<fftwl_complex*>(y), NULL, strideY, distY,
                                      direction, FFTW_PLANNER_FLAG);
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
dft_multiple(IndexType n, IndexType m,
             std::complex<__float128> *x, IndexType strideX, IndexType distX,
             std::complex<__float128> *y, IndexType strideY, IndexType distY,
             DFTDirection direction)
{
    CXXBLAS_DEBUG_OUT("dft_multiple [FFTW interface, quad]");

#   ifdef FFTW_WISDOM_IMPORT
        ASSERT(strcmp(FFTW_WISDOM_FILENAME,""));
        fftwq_import_wisdom_from_filename(FFTW_WISDOM_FILENAME);
#   endif

    fftwq_plan p = fftwq_plan_many_dft(1, &n, m,
                                      reinterpret_cast<fftwq_complex*>(x), NULL, strideX, distX,
                                      reinterpret_cast<fftwq_complex*>(y), NULL, strideY, distY,
                                      direction, FFTW_PLANNER_FLAG);
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

#endif // PLAYGROUND_CXXDFT_MULTIPLE_TCC
