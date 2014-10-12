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

#ifndef FLENS_STORAGE_BANDSTORAGE_HASBANDSTORAGE_H
#define FLENS_STORAGE_BANDSTORAGE_HASBANDSTORAGE_H 1

#include <flens/auxiliary/hasengine.h>
#include <flens/storage/bandstorage/isbandstorage.h>

namespace flens {

template <typename T, bool>
struct GetBandEngine_
{
    typedef void Type;
};

template <typename T>
struct GetBandEngine_<T, true>
{
    typedef typename T::Engine Type;
};

template <typename T>
struct HasBandStorage
{
    typedef typename GetBandEngine_<T, HasEngine<T>::value>::Type  Engine;
    static const bool value = IsBandStorage<Engine>::value;
};

} // namespace flens

#endif // FLENS_STORAGE_BANDSTORAGE_HASBANDSTORAGE_H
