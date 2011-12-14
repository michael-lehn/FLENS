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

#ifndef FLENS_LAPACK_AUX_CONVERT_H
#define FLENS_LAPACK_AUX_CONVERT_H 1

#include <flens/typedefs.h>
#include <flens/lapack/typedefs.h>
#include <flens/lapack/aux/lascl.h>
#include <flens/lapack/eig/bak.h>
#include <flens/lapack/eig/bal.h>
#include <flens/lapack/eig/evx.h>
#include <flens/lapack/eig/hseqr.h>
#include <flens/lapack/eig/trevc.h>
#include <flens/lapack/eig/trsna.h>

namespace flens { namespace lapack {

//-- convert: enum Transpose
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,Transpose>::value, char>::Type
    getF77LapackChar(ENUM enumValue);

//-- convert: enum Side
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,Side>::value, char>::Type
    getF77LapackChar(ENUM enumValue);

//-- convert: enum Norm
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,Norm>::value, char>::Type
    getF77LapackChar(ENUM enumValue);

//-- convert: enum BALANCE::Balance
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,BALANCE::Balance>::value, char>::Type
    getF77LapackChar(ENUM enumValue);

//-- convert: enum SENSE::Sense
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,SENSE::Sense>::value, char>::Type
    getF77LapackChar(ENUM enumValue);

//-- convert: enum HSEQR::Job
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,HSEQR::Job>::value, char>::Type
    getF77LapackChar(ENUM enumValue);

//-- convert: enum HSEQR::ComputeZ
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,HSEQR::ComputeZ>::value, char>::Type
    getF77LapackChar(ENUM enumValue);

//-- convert: enum TREVC::Job
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,TREVC::Job>::value, char>::Type
    getF77LapackChar(ENUM enumValue);

//-- convert: enum TRSNA::Job
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,TRSNA::Job>::value, char>::Type
    getF77LapackChar(ENUM enumValue);

//-- convert: enum TRSNA::HowMany
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,TRSNA::HowMany>::value, char>::Type
    getF77LapackChar(ENUM enumValue);

//-- convert: enum LASCL::Type
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,LASCL::Type>::value, char>::Type
    getF77LapackChar(ENUM enumValue);

//------------------------------------------------------------------------------

//-- convert: char to enum BALANCE::Balance
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,BALANCE::Balance>::value,
                        BALANCE::Balance>::Type
    getFlensLapackEnum(char trans);

//-- convert: char to enum SENSE::Sense
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,SENSE::Sense>::value,
                        SENSE::Sense>::Type
    getFlensLapackEnum(char trans);

//-- convert: char to enum Transpose
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,Transpose>::value, Transpose>::Type
    getFlensLapackEnum(char trans);

//-- convert: char to enum Side
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,Side>::value, Side>::Type
    getFlensLapackEnum(char side);

//-- convert: char to enum Norm
template <typename ENUM>
    typename RestrictTo<IsSame<ENUM,Norm>::value, Norm>::Type
    getFlensLapackEnum(char norm);

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_CONVERT_H
