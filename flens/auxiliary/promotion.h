/*
 *   Copyright (c) 2010, Michael Lehn
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

#ifndef FLENS_AUXILIARY_PROMOTION_H
#define FLENS_AUXILIARY_PROMOTION_H 1

#include <flens/auxiliary/ifelseif.h>
#include <flens/auxiliary/isconvertible.h>

namespace flens {

template <typename... Args>
struct Promotion
{
};

template <typename T>
struct Promotion<T>
{
    typedef T Type;
};

template <typename T>
struct Promotion<T, T>
{
    typedef T Type;
};

template <typename T, typename... Args>
struct Promotion<T, Args...>
{
    typedef typename Promotion<Args... >::Type  ArgsT;
    typedef typename Promotion<T, ArgsT>::Type  Type;
};

template <typename T1, typename T2>
struct Promotion<T1, T2>
{
    typedef typename IfElseIf<sizeof(T2)<=sizeof(T1), T1,
                              sizeof(T1)<=sizeof(T2), T2>::Type     Type;
};

//-- int, float ----------------------------------------------------------------
template <>
struct Promotion<int, float>
{
    typedef float  Type;
};

template <>
struct Promotion<float, int>
{
    typedef float  Type;
};

//-- long, float ---------------------------------------------------------------
template <>
struct Promotion<long, float>
{
    typedef float  Type;
};

template <>
struct Promotion<float, long>
{
    typedef float  Type;
};

//-- int, double ---------------------------------------------------------------
template <>
struct Promotion<int, double>
{
    typedef double  Type;
};

template <>
struct Promotion<double, int>
{
    typedef double  Type;
};

//-- long, double --------------------------------------------------------------
template <>
struct Promotion<long, double>
{
    typedef double  Type;
};

template <>
struct Promotion<double, long>
{
    typedef double  Type;
};

} // namespace flens

#endif // FLENS_AUXILIARY_PROMOTION_H
