/*
 *   Copyright (c) 2009, Michael Lehn
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

#ifndef FLENS_BLAS_LEVEL1_SCAL_H
#define FLENS_BLAS_LEVEL1_SCAL_H 1

#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//-- scal
template <typename ALPHA, typename VY>
    typename RestrictTo<IsDenseVector<VY>::value,
             void>::Type
    scal(const ALPHA &alpha, VY &&y);

//-- gescal
template <typename ALPHA, typename MB>
    typename RestrictTo<IsGeMatrix<MB>::value,
             void>::Type
    scal(const ALPHA &alpha, MB &&B);

//-- gbscal
template <typename ALPHA, typename MB>
    typename RestrictTo<IsGbMatrix<MB>::value,
             void>::Type
    scal(const ALPHA &alpha, MB &&B);

//-- hbscal
template <typename ALPHA, typename MB>
    typename RestrictTo<IsHbMatrix<MB>::value,
             void>::Type
    scal(const ALPHA &alpha, MB &&B);

//-- hpscal
template <typename ALPHA, typename MB>
    typename RestrictTo<IsHpMatrix<MB>::value,
             void>::Type
    scal(const ALPHA &alpha, MB &&B);

//-- sbscal
template <typename ALPHA, typename MB>
    typename RestrictTo<IsSbMatrix<MB>::value,
             void>::Type
    scal(const ALPHA &alpha, MB &&B);

//-- spscal
template <typename ALPHA, typename MB>
    typename RestrictTo<IsSpMatrix<MB>::value,
             void>::Type
    scal(const ALPHA &alpha, MB &&B);

//-- tbscal
template <typename ALPHA, typename MB>
    typename RestrictTo<IsTbMatrix<MB>::value,
             void>::Type
    scal(const ALPHA &alpha, MB &&B);

//-- tpscal
template <typename ALPHA, typename MB>
    typename RestrictTo<IsTpMatrix<MB>::value,
             void>::Type
    scal(const ALPHA &alpha, MB &&B);

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_SCAL_H
