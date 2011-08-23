/*
 *   Copyright (c) 2011, Michael Lehn
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

#ifndef FLENS_STORAGE_FULLSTORAGE_TRAPEZOIDALFILL_TCC
#define FLENS_STORAGE_FULLSTORAGE_TRAPEZOIDALFILL_TCC 1

#include <cxxblas/typedefs.h>

namespace flens {

template <typename IndexType, typename T>
void
trapezoidalFill(cxxblas::StorageOrder order, cxxblas::StorageUpLo upLo,
                const T &value,
                IndexType m, IndexType n, T *data, IndexType ld)
{
    // fill the strictly upper trapezoidal part
    if (upLo==cxxblas::Upper) {
        if (order==cxxblas::RowMajor) {
            for (IndexType i=0; i<min(m,n); ++i, data+=ld+1) {
                std::fill_n(data+1, n-i-1, value);
            }
        }
        if (order==cxxblas::ColMajor) {
            ++data;
            for (IndexType j=1; j<n; ++j, data+=ld) {
                std::fill_n(data, min(j,m), value);
            }
        }
    }
    // fill the strictly lower trapezoidal part
    if (upLo==cxxblas::Lower) {
        if (order==cxxblas::RowMajor) {
            data += ld;
            for (IndexType i=1; i<m; ++i, data+=ld) {
                std::fill_n(data, min(i, n), value);
            }
        }
        if (order==cxxblas::ColMajor) {
            for (IndexType j=0; j<min(m, n); ++j, data+=ld+1) {
                std::fill_n(data+1, m-j-1, value);
            }
        }
    }
}

} // namespace flens

#endif // FLENS_STORAGE_FULLSTORAGE_TRAPEZOIDALFILL_TCC
