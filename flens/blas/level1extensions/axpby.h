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

#ifndef FLENS_BLAS_LEVEL1_AXPBY_H
#define FLENS_BLAS_LEVEL1_AXPBY_H 1

#include <cxxblas/cxxblas.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/storage/storage.h>
#include <flens/typedefs.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//-- BLAS Level 1 extensions ---------------------------------------------------

//-- axpby
template <typename ALPHA, typename VX, typename BETA, typename VY>
    typename RestrictTo<IsDenseVector<VX>::value
                     && IsDenseVector<VY>::value,
             void>::Type
    axpby(const ALPHA &alpha, const VX &x, const BETA &beta, VY &&y);


//-- axpbby
template <typename ALPHA, typename VX, typename BETA, typename VY>
    typename RestrictTo<IsTinyVector<VX>::value
                     && IsTinyVector<VY>::value,
             void>::Type
    axpby(const ALPHA &alpha, const VX &x, const BETA &beta, VY &&y);

//== GeneralMatrix

//-- diagaxpby
template <typename ALPHA, typename MA, typename BETA, typename MB>
    typename RestrictTo<IsDiagMatrix<MA>::value
                     && IsDiagMatrix<MB>::value,
    void>::Type
    axpby(Transpose trans, const ALPHA &alpha, const MA &A,
          const BETA &beta, MB &&B);

//-- gbaxpby
template <typename ALPHA, typename MA, typename BETA, typename MB>
    typename RestrictTo<IsGbMatrix<MA>::value
                     && IsGbMatrix<MB>::value,
    void>::Type
    axpby(Transpose trans, const ALPHA &alpha, const MA &A,
          const BETA &beta, MB &&B);

//-- geaxpby
template <typename ALPHA, typename MA, typename BETA, typename MB>
    typename RestrictTo<IsGeMatrix<MA>::value
                     && IsGeMatrix<MB>::value,
             void>::Type
    axpby(Transpose trans, const ALPHA &alpha, const MA &A,
          const BETA &beta, MB &&B);

//-- (tiny) geaxpby
template <typename ALPHA, typename MA, typename BETA, typename MB>
    typename RestrictTo<IsGeTinyMatrix<MA>::value
                     && IsGeTinyMatrix<MB>::value,
             void>::Type
    axpby(Transpose trans, const ALPHA &alpha, const MA &A,
          const BETA &beta, MB &&B);

//== HermitianMatrix

//-- hbaxpby
template <typename ALPHA, typename MA, typename BETA, typename MB>
    typename RestrictTo<IsHbMatrix<MA>::value
                    && IsHbMatrix<MB>::value,
             void>::Type
    axpby(Transpose trans, const ALPHA &alpha, const MA &A,
          const BETA &beta, MB &&B);

//-- heaxpby
template <typename ALPHA, typename MA, typename BETA, typename MB>
    typename RestrictTo<IsHeMatrix<MA>::value
                     && IsHeMatrix<MB>::value,
             void>::Type
    axpby(Transpose trans, const ALPHA &alpha, const MA &A,
          const BETA &beta, MB &&B);

//-- hpaxpby
template <typename ALPHA, typename MA, typename BETA, typename MB>
    typename RestrictTo<IsHpMatrix<MA>::value
                     && IsHpMatrix<MB>::value,
             void>::Type
    axpby(Transpose trans, const ALPHA &alpha, const MA &A,
          const BETA &beta, MB &&B);


//== SymmetricMatrix

//-- sbaxpby
template <typename ALPHA, typename MA, typename BETA, typename MB>
    typename RestrictTo<IsSbMatrix<MA>::value
                     && IsSbMatrix<MB>::value,
             void>::Type
    axpby(Transpose trans, const ALPHA &alpha, const MA &A,
          const BETA &beta, MB &&B);

//-- spaxpby
template <typename ALPHA, typename MA, typename BETA, typename MB>
    typename RestrictTo<IsSpMatrix<MA>::value
                     && IsSpMatrix<MB>::value,
             void>::Type
    axpby(Transpose trans, const ALPHA &alpha, const MA &A,
          const BETA &beta, MB &&B);

//-- syaxpby
template <typename ALPHA, typename MA, typename BETA, typename MB>
    typename RestrictTo<IsSyMatrix<MA>::value
                     && IsSyMatrix<MB>::value,
             void>::Type
    axpby(Transpose trans, const ALPHA &alpha, const MA &A,
          const BETA &beta, MB &&B);


//== TriangularMatrix

//-- tbaxpby
template <typename ALPHA, typename MA, typename BETA, typename MB>
    typename RestrictTo<IsTbMatrix<MA>::value
                     && IsTbMatrix<MB>::value,
             void>::Type
    axpby(Transpose trans, const ALPHA &alpha, const MA &A,
          const BETA &beta, MB &&B);

//-- traxpby
template <typename ALPHA, typename MA, typename BETA, typename MB>
    typename RestrictTo<IsTrMatrix<MA>::value
                     && IsTrMatrix<MB>::value,
             void>::Type
    axpby(Transpose trans, const ALPHA &alpha, const MA &A,
          const BETA &beta, MB &&B);

//-- tpaxpby
template <typename ALPHA, typename MA, typename BETA, typename MB>
    typename RestrictTo<IsTpMatrix<MA>::value
                     && IsTpMatrix<MB>::value,
             void>::Type
    axpby(Transpose trans, const ALPHA &alpha, const MA &A,
          const BETA &beta, MB &&B);

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_AXPBY_H
