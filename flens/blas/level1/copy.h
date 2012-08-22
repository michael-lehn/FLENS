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

#ifndef FLENS_BLAS_LEVEL1_COPY_H
#define FLENS_BLAS_LEVEL1_COPY_H 1

#include <cxxblas/cxxblas.h>
#include <flens/auxiliary/auxiliary.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//-- BLAS Level 1 --------------------------------------------------------------

//-- copy
template <typename VX, typename VY>
    typename RestrictTo<IsDenseVector<VX>::value
                     && IsDenseVector<VY>::value,
             void>::Type
    copy(const VX &x, VY &&y);

//-- BLAS Level 1 extensions ---------------------------------------------------

//-- gbcopy
template <typename MA, typename MB>
    typename RestrictTo<IsGbMatrix<MA>::value
                     && IsGeMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);

//-- gecopy
template <typename MA, typename MB>
    typename RestrictTo<IsGeMatrix<MA>::value
                     && IsGeMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);

//-- hbcopy
template <typename MA, typename MB>
    typename RestrictTo<IsHbMatrix<MA>::value
                     && IsHbMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);

//-- hpcopy
template <typename MA, typename MB>
    typename RestrictTo<IsHpMatrix<MA>::value
                     && IsHpMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);
    
//-- tbcopy
template <typename MA, typename MB>
    typename RestrictTo<IsTbMatrix<MA>::value
                     && IsTbMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);

//-- trcopy
template <typename MA, typename MB>
    typename RestrictTo<IsTrMatrix<MA>::value
                     && IsTrMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);

//-- tpcopy
template <typename MA, typename MB>
    typename RestrictTo<IsTpMatrix<MA>::value
                     && IsTpMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);

//-- sbcopy
template <typename MA, typename MB>
    typename RestrictTo<IsSbMatrix<MA>::value
                     && IsSbMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);

//-- spcopy
template <typename MA, typename MB>
    typename RestrictTo<IsSpMatrix<MA>::value
                     && IsSpMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);

//-- sycopy
template <typename MA, typename MB>
    typename RestrictTo<IsSyMatrix<MA>::value
                     && IsSyMatrix<MB>::value,
             void>::Type
    copy(const MA &A, MB &&B);

//-- Sparse BLAS extensions ----------------------------------------------------
    
//-- copy: GeCoordMatrix -> GeCCSMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsGeCoordMatrix<MA>::value
                     && IsGeCCSMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);

//-- copy: GeCoordMatrix -> GeCRSMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsGeCoordMatrix<MA>::value
                     && IsGeCRSMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);

//-- copy: SyCoordMatrix -> SyCCSMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsSyCoordMatrix<MA>::value
                     && IsSyCCSMatrix<MB>::value,
             void>::Type
    copy(const MA &A, MB &&B);

//-- copy: SyCoordMatrix -> SyCRSMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsSyCoordMatrix<MA>::value
                     && IsSyCRSMatrix<MB>::value,
             void>::Type
    copy(const MA &A, MB &&B);

//-- convenience extensions ----------------------------------------------------
    
//-- copy: HbMatrix -> GbMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsHbMatrix<MA>::value
                     && IsGbMatrix<MB>::value,
             void>::Type
    copy(const MA &A, MB &&B);
    
//-- copy: TbMatrix -> GbMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsTbMatrix<MA>::value
                     && IsGbMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);
    
//-- copy: TrMatrix -> GeMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsTrMatrix<MA>::value
                     && IsGeMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);
    
//-- copy: SbMatrix -> GbMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsSbMatrix<MA>::value
                     && IsGbMatrix<MB>::value,
             void>::Type
    copy(const MA &A, MB &&B);
    
//-- copy: SyMatrix -> GeMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsSyMatrix<MA>::value
                     && IsGeMatrix<MB>::value,
             void>::Type
    copy(const MA &A, MB &&B);

//-- copy: GeCoordMatrix -> GeMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsGeCoordMatrix<MA>::value
                     && IsGeMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);

//-- copy: GeCCSMatrix -> GeMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsGeCCSMatrix<MA>::value
                     && IsGeMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);

//-- copy: GeCRSMatrix -> GeMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsGeCRSMatrix<MA>::value
                     && IsGeMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const MA &A, MB &&B);

//-- copy: SyCoordMatrix -> SyMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsSyCoordMatrix<MA>::value
                     && IsSyMatrix<MB>::value,
             void>::Type
    copy(const MA &A, MB &&B);

//-- copy: SyCCSMatrix -> SyMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsSyCCSMatrix<MA>::value
                     && IsSyMatrix<MB>::value,
             void>::Type
    copy(const MA &A, MB &&B);

//-- copy: SyCRSMatrix -> SyMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsSyCRSMatrix<MA>::value
                     && IsSyMatrix<MB>::value,
             void>::Type
    copy(const MA &A, MB &&B);

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_COPY_H
