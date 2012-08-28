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
#include <flens/auxiliary/auxiliary.h>
#include <flens/typedefs.h>

namespace flens {

template <typename IndexType, typename T>
void
trapezoidalFill(StorageOrder order, StorageUpLo upLo,
                const T &value,
                IndexType m, IndexType n, T *data, IndexType ld)
{
    using std::min;

    // fill the upper trapezoidal part
    if (upLo==Upper) {
        if (order==RowMajor) {
            for (IndexType i=0; i<min(m,n); ++i, data+=ld+1) {
                std::fill_n(data, n-i, value);
            }
        }
        if (order==ColMajor) {
            for (IndexType j=0; j<n; ++j, data+=ld) {
                std::fill_n(data, min(j+1,m), value);
            }
        }
    }
    // fill the lower trapezoidal part
    if (upLo==Lower) {
        if (order==RowMajor) {
            for (IndexType i=0; i<m; ++i, data+=ld) {
                std::fill_n(data, min(i+1, n), value);
            }
        }
        if (order==ColMajor) {
            for (IndexType j=0; j<min(m, n); ++j, data+=ld+1) {
                std::fill_n(data, m-j, value);
            }
        }
    }
}

template <typename IndexType, typename T>
void
trapezoidalFillRandom(StorageOrder order, StorageUpLo upLo,
                      IndexType m, IndexType n, T *data, IndexType ld)
{
    using std::min;

    // fill the upper trapezoidal part
    if (upLo==Upper) {
        if (order==RowMajor) {
            for (IndexType i=0; i<min(m,n); ++i, data+=ld+1) {
                for (IndexType j=0; j<n-i;++j) {
                    data[j] = randomValue<T>();
                }
            }
        }
        if (order==ColMajor) {
            for (IndexType j=0; j<n; ++j, data+=ld) {
                for (IndexType i=0; i<min(j+1,m); ++i) {
                    data[i] = randomValue<T>();
                }
            }
        }
    }
    // fill the lower trapezoidal part
    if (upLo==Lower) {
        if (order==RowMajor) {
            for (IndexType i=0; i<m; ++i, data+=ld) {
                for (IndexType j=0; j<min(i+1,n); ++j) {
                    data[j] = randomValue<T>();
                }
            }
        }
        if (order==ColMajor) {
            for (IndexType j=0; j<min(m, n); ++j, data+=ld+1) {
                for (IndexType i=0; i<m-j; ++i) {
                    data[i] = randomValue<T>();
                }
            }
        }
    }
}

} // namespace flens

#endif // FLENS_STORAGE_FULLSTORAGE_TRAPEZOIDALFILL_TCC
