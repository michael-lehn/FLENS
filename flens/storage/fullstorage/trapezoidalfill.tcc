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

#include <flens/auxiliary/auxiliary.h>
#include <flens/typedefs.h>

namespace flens {

template <typename IndexType, typename T>
void
trapezoidalFill(IndexType m, IndexType n, const T &value,
                StorageUpLo upLo,
                T *data, IndexType incRow, IndexType incCol)
{
    if (upLo==Upper) {
        trapezoidalFillRandom(n, m, Lower, data, incCol, incRow);
        return;
    }

    if (incRow<incCol) {
        const IndexType k = std::min(m, n);
        for (IndexType j=0; j<k; ++j) {
            for (IndexType i=j; i<m; ++i) {
                data[i*incRow+j*incCol] = value;
            }
        }
    } else {
        for (IndexType i=0; i<m; ++i) {
            const IndexType k = std::min(i+1, n);
            for (IndexType j=0; j<k; ++j) {
                data[i*incRow+j*incCol] = value;
            }
        }
    }
}

template <typename IndexType, typename T>
void
trapezoidalFillRandom(IndexType m, IndexType n,
                      StorageUpLo upLo,
                      T *data, IndexType incRow, IndexType incCol)
{
    if (upLo==Upper) {
        trapezoidalFillRandom(n, m, Lower, data, incCol, incRow);
        return;
    }

    if (incRow<incCol) {
        const IndexType k = std::min(m, n);
        for (IndexType j=0; j<k; ++j) {
            for (IndexType i=j; i<m; ++i) {
                data[i*incRow+j*incCol] = randomValue<T>();
            }
        }
    } else {
        for (IndexType i=0; i<m; ++i) {
            const IndexType k = std::min(i+1, n);
            for (IndexType j=0; j<k; ++j) {
                data[i*incRow+j*incCol] = randomValue<T>();
            }
        }
    }
}

} // namespace flens

#endif // FLENS_STORAGE_FULLSTORAGE_TRAPEZOIDALFILL_TCC
