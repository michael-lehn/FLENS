/*
 *   Copyright (c) 2012, Michael Lehn
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

#ifndef FLENS_STORAGE_TINYARRAY_TINYARRAY_TCC
#define FLENS_STORAGE_TINYARRAY_TINYARRAY_TCC 1

#include <memory>
#include <cxxblas/cxxblas.h>
#include <flens/storage/tinyarray/tinyarray.h>

namespace flens {

template <typename T, int n, int ib>
TinyArray<T,n,ib>::TinyArray()
{
}

template <typename T, int n, int ib>
TinyArray<T,n,ib>::TinyArray(const TinyArray &rhs)
{
    cxxblas::copy<n, T, T>(rhs.data(), data());
}

template <typename T, int n, int ib>
TinyArray<T,n,ib>::~TinyArray()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, int n, int ib>
const typename TinyArray<T,n,ib>::ElementType &
TinyArray<T,n,ib>::operator()(IndexType index) const
{
    const T *data = reinterpret_cast<const T *>(_data) - firstIndex;
    return data[index];
}

template <typename T, int n, int ib>
typename TinyArray<T,n,ib>::ElementType &
TinyArray<T,n,ib>::operator()(IndexType index)
{
    T *data = reinterpret_cast<T *>(_data) - firstIndex;
    return data[index];
}

//-- methods -------------------------------------------------------------------

template <typename T, int n, int ib>
const typename TinyArray<T,n,ib>::ElementType *
TinyArray<T,n,ib>::data() const
{
    return reinterpret_cast<const T *>(_data);
}

template <typename T, int n, int ib>
typename TinyArray<T,n,ib>::ElementType *
TinyArray<T,n,ib>::data()
{
    return reinterpret_cast<T *>(_data);
}

template <typename T, int n, int ib>
void
TinyArray<T,n,ib>::fill(const ElementType &value)
{
    std::fill_n(data(), length, value);
}

} // namespace flens

#endif // FLENS_STORAGE_TINYARRAY_TINYARRAY_TCC
