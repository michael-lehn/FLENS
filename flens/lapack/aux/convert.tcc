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

#ifndef FLENS_LAPACK_AUX_CONVERT_TCC
#define FLENS_LAPACK_AUX_CONVERT_TCC 1

#include <flens/aux/aux.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//-- convert: enum Transpose
template <typename ENUM>
typename RestrictTo<IsSame<ENUM,Transpose>::value, char>::Type
getF77LapackChar(ENUM enumValue)
{
    if (enumValue==NoTrans) {
        return 'N';
    } else if (enumValue==Trans) {
        return 'T';
    } else if (enumValue==Conj) {
        return 'R';
    } else if (enumValue==ConjTrans) {
        return 'H';
    } else {
        ASSERT(0);
        return '?';
    }
}

//-- convert: enum Side
template <typename ENUM>
typename RestrictTo<IsSame<ENUM,Side>::value, char>::Type
getF77LapackChar(ENUM enumValue)
{
    if (enumValue==Left) {
        return 'L';
    } else if (enumValue==Right) {
        return 'R';
    } else {
        ASSERT(0);
        return '?';
    }
}

//-- convert: enum Norm
template <typename ENUM>
typename RestrictTo<IsSame<ENUM,Norm>::value, char>::Type
getF77LapackChar(ENUM enumValue)
{
    if (enumValue==MaximumNorm) {
        return 'M';
    } else if (enumValue==OneNorm) {
        return 'O';
    } else if (enumValue==InfinityNorm) {
        return 'I';
    } else if (enumValue==FrobeniusNorm) {
        return 'F';
    } else {
        ASSERT(0);
    }
}

//-- convert: enum BALANCE::Balance
template <typename ENUM>
typename RestrictTo<IsSame<ENUM,BALANCE::Balance>::value, char>::Type
getF77LapackChar(ENUM enumValue)
{
    if (enumValue==BALANCE::None) {
        return 'N';
    } else if (enumValue==BALANCE::PermuteOnly) {
        return 'P';
    } else if (enumValue==BALANCE::ScaleOnly) {
        return 'S';
    } else if (enumValue==BALANCE::Both) {
        return 'B';
    } else {
        ASSERT(0);
    }
}

//-- convert: enum SENSE::Sense
template <typename ENUM>
typename RestrictTo<IsSame<ENUM,SENSE::Sense>::value, char>::Type
getF77LapackChar(ENUM enumValue)
{
    if (enumValue==SENSE::None) {
        return 'N';
    } else if (enumValue==SENSE::EigenvaluesOnly) {
        return 'E';
    } else if (enumValue==SENSE::InvariantSubspaceOnly) {
        return 'V';
    } else if (enumValue==SENSE::Both) {
        return 'B';
    } else {
        ASSERT(0);
    }
}

//-- convert: enum HSEQR::Job
template <typename ENUM>
typename RestrictTo<IsSame<ENUM,HSEQR::Job>::value, char>::Type
getF77LapackChar(ENUM enumValue)
{
    if (enumValue==HSEQR::Eigenvalues) {
        return 'E';
    } else if (enumValue==HSEQR::Schur) {
        return 'S';
    } else {
        ASSERT(0);
    }
}

//-- convert: enum HSEQR::ComputeZ
template <typename ENUM>
typename RestrictTo<IsSame<ENUM,HSEQR::ComputeZ>::value, char>::Type
getF77LapackChar(ENUM enumValue)
{
    if (enumValue==HSEQR::No) {
        return 'N';
    } else if (enumValue==HSEQR::Init) {
        return 'I';
    } else if (enumValue==HSEQR::NoInit) {
        return 'V';
    } else {
        ASSERT(0);
    }
}

//-- convert: enum TREVC::Job
template <typename ENUM>
typename RestrictTo<IsSame<ENUM,TREVC::Job>::value, char>::Type
getF77LapackChar(ENUM enumValue)
{
    if (enumValue==TREVC::All) {
        return 'A';
    } else if (enumValue==TREVC::Backtransform) {
        return 'B';
    } else if (enumValue==TREVC::Selected) {
        return 'S';
    } else {
        ASSERT(0);
    }
}

//-- convert: enum TRSNA::Job
template <typename ENUM>
typename RestrictTo<IsSame<ENUM,TRSNA::Job>::value, char>::Type
getF77LapackChar(ENUM enumValue)
{
    if (enumValue==TRSNA::EigenvaluesOnly) {
        return 'E';
    } else if (enumValue==TRSNA::EigenvectorsOnly) {
        return 'V';
    } else if (enumValue==TRSNA::Both) {
        return 'B';
    } else {
        ASSERT(0);
    }
}

//-- convert: enum TRSNA::HowMany
template <typename ENUM>
typename RestrictTo<IsSame<ENUM,TRSNA::HowMany>::value, char>::Type
getF77LapackChar(ENUM enumValue)
{
    if (enumValue==TRSNA::All) {
        return 'A';
    } else if (enumValue==TRSNA::Selected) {
        return 'S';
    } else {
        ASSERT(0);
    }
}

//-- convert: enum LASCL::Type
template <typename ENUM>
typename RestrictTo<IsSame<ENUM,LASCL::Type>::value, char>::Type
getF77LapackChar(ENUM enumValue)
{
    if (enumValue==LASCL::FullMatrix) {
        return 'G';
    } else if (enumValue==LASCL::LowerTriangular) {
        return 'L';
    } else if (enumValue==LASCL::UpperTriangular) {
        return 'U';
    } else if (enumValue==LASCL::UpperHessenberg) {
        return 'H';
    } else if (enumValue==LASCL::SymmetricLowerBand) {
        return 'B';
    } else if (enumValue==LASCL::SymmetricUpperBand) {
        return 'Q';
    } else if (enumValue==LASCL::GeneralBand) {
        return 'Z';
    } else {
        ASSERT(0);
    }
}

//------------------------------------------------------------------------------
//-- convert: char to enum BALANCE::Balance
template <typename ENUM>
typename RestrictTo<IsSame<ENUM,BALANCE::Balance>::value,
                    BALANCE::Balance>::Type
getFlensLapackEnum(char balance)
{
    if (balance=='N') {
        return BALANCE::None;
    } else if (balance=='P') {
        return BALANCE::PermuteOnly;
    } else if (balance=='S') {
        return BALANCE::ScaleOnly;
    } else if (balance=='B') {
        return BALANCE::Both;
    } else {
        ASSERT(0);
    }
}

//-- convert: char to enum SENSE::Sense
template <typename ENUM>
typename RestrictTo<IsSame<ENUM,SENSE::Sense>::value, SENSE::Sense>::Type
getFlensLapackEnum(char sense)
{
    if (sense=='N') {
        return SENSE::None;
    } else if (sense=='E') {
        return SENSE::EigenvaluesOnly;
    } else if (sense=='V') {
        return SENSE::InvariantSubspaceOnly;
    } else if (sense=='B') {
        return SENSE::Both;
    } else {
        ASSERT(0);
        // this should not happen ...
        return SENSE::None;
    }
}

//-- convert: char to enum Transpose
template <typename ENUM>
typename RestrictTo<IsSame<ENUM,Transpose>::value, Transpose>::Type
getFlensLapackEnum(char trans)
{
    if (trans=='N') {
        return NoTrans;
    } else if (trans=='T') {
        return Trans;
    } else if (trans=='C') {
        return ConjTrans;
    } else if (trans=='R') {
        return Conj;
    }
    ASSERT(0);
}

//-- convert: char to enum Side
template <typename ENUM>
typename RestrictTo<IsSame<ENUM,Side>::value, Side>::Type
getFlensLapackEnum(char side)
{
    if (side=='L') {
        return Left;
    } else if (side=='R') {
        return Right;
    }
    ASSERT(0);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_CONVERT_TCC
