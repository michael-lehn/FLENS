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

#ifndef FLENS_AUXILIARY_MAX_TCC
#define FLENS_AUXILIARY_MAX_TCC 1

#include <algorithm>
#include <flens/auxiliary/max.h>
#include <flens/auxiliary/restrictto.h>

namespace flens {

//
//  Variant 1:
//  All arguments have same type.  Function returns a const reference.
//
template <typename T>
const T &
max(const T &a)
{
    return a;
}

template <typename T>
const T &
max(const T &a, const T &b)
{
    return std::max(a, b);
}

template <typename T, typename ...Args>
const typename RestrictTo<IsSame<T,Args...>::value, T>::Type &
max(const T &a, const T &b, const Args &...args)
{
    return flens::max(flens::max(a, b), args...);
}

//
//  Variant 2:
//  Arguments have different types.  Function returns a copy
//
template <typename T1, typename T2>
const typename Promotion<T1, T2>::Type
max(const T1 &a, const T2 &b)
{
    return (a>b) ? a : b;
}

template <typename T1, typename T2, typename ...Args>
const typename RestrictTo<!IsSame<T1, T2, Args...>::value,
                          typename Promotion<T1, T2, Args...>::Type>::Type
max(const T1 &a, const T2 &b, const Args & ...args)
{
    return max(max(a, b), args...);
}

} // namespace flens

#endif // FLENS_AUXILIARY_MAX_TCC
