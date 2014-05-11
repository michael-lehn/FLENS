/*
 *   Copyright (c) 2014, Michael Lehn
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

/* Based on
 *
       SUBROUTINE ZSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_STEQR_H
#define FLENS_LAPACK_IMPL_STEQR_H 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/lapack/typedefs.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {


//== steqr =====================================================================

namespace STEQR {

    enum ComputeZ {
        No     = 'N',   // Compute eigenvalues only.
        Orig   = 'V',   // Compute eigenvalues and eigenvectors of the original
                        // Hermitian matrix.  On entry, Z must contain the
                        // unitary matrix used to reduce the original matrix
                        // to tridiagonal form.
        Tri    = 'I',   // Compute eigenvalues and eigenvectors of the
                        // tridiagonal matrix.  Z is initialized to the identity
                        // matrix.
    };

}

template <typename VD, typename VE, typename MZ, typename VWORK>
    typename RestrictTo<IsRealDenseVector<VD>::value
                     && IsRealDenseVector<VE>::value
                     && IsComplexGeMatrix<MZ>::value
                     && IsRealDenseVector<VWORK>::value,
             typename RemoveRef<VD>::Type::IndexType>::Type
    steqr(STEQR::ComputeZ  compZ,
          VD               &&d,
          VE               &&e,
          MZ               &&Z,
          VWORK            &&work);

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_STEQR_H
