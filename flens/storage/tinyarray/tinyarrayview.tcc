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

#ifndef FLENS_STORAGE_TINYARRAY_TINYARRAYVIEW_TCC
#define FLENS_STORAGE_TINYARRAY_TINYARRAYVIEW_TCC 1

#include <flens/storage/tinyarray/tinyarrayview.h>

namespace flens {

template <typename T, int n, int inc, int indexBase>
TinyArrayView<T,n,inc,indexBase>::TinyArrayView(ElementType *data)
    : data_(data)
{
}

template <typename T, int n, int inc, int indexBase>
TinyArrayView<T,n,inc,indexBase>::~TinyArrayView()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, int n, int inc, int indexBase>
const typename TinyArrayView<T,n,inc,indexBase>::ElementType &
TinyArrayView<T,n,inc,indexBase>::operator()(IndexType index) const
{
    return data_[inc*(index-indexBase)];
}

template <typename T, int n, int inc, int indexBase>
typename TinyArrayView<T,n,inc,indexBase>::ElementType &
TinyArrayView<T,n,inc,indexBase>::operator()(IndexType index)
{
    return data_[inc*(index-indexBase)];
}

//-- methods -------------------------------------------------------------------

template <typename T, int n, int inc, int indexBase>
const typename TinyArrayView<T,n,inc,indexBase>::ElementType *
TinyArrayView<T,n,inc,indexBase>::data() const
{
    return data_;
}

template <typename T, int n, int inc, int indexBase>
typename TinyArrayView<T,n,inc,indexBase>::ElementType *
TinyArrayView<T,n,inc,indexBase>::data()
{
    return data_;
}

template <typename T, int n, int inc, int indexBase>
void
TinyArrayView<T,n,inc,indexBase>::fill(const ElementType &value)
{
    for (int i=0, I=0; i<n; ++i, I+=inc) {
        data_[I] = value;
    }
}

} // namespace flens

#endif // FLENS_STORAGE_TINYARRAY_TINYARRAYVIEW_TCC
