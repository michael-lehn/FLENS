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

/* Based on
 *
       SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
 *
 *  -- LAPACK auxiliary routine (version 3.3.0) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2010
 */

#ifndef FLENS_LAPACK_IMPL_LASCL_H
#define FLENS_LAPACK_IMPL_LASCL_H 1


#include <flens/auxiliary/auxiliary.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

namespace LASCL {

    enum Type {
        FullMatrix = 'G',         // = 'G':  A is a full matrix.
        LowerTriangular = 'L',    // = 'L':  A is a lower triangular matrix.
        UpperTriangular = 'U',    // = 'U':  A is an upper triangular matrix.
        UpperHessenberg = 'H',    // = 'H':  A is an upper Hessenberg matrix.
        SymmetricLowerBand = 'B', // = 'B':  A is a symmetric band matrix with
                                  //         lower bandwidth KL and upper
                                  //         bandwidth KU and with the only
                                  //         the lower half stored.
        SymmetricUpperBand = 'Q', // = 'Q':  A is a symmetric band matrix with
                                  //         lower bandwidth KL and upper
                                  //         bandwidth KU and with the only the
                                  //         upper half stored.
        GeneralBand = 'Z'         // = 'Z':  A is a band matrix with lower
                                  //         bandwidth KL and upper
                                  //         bandwidth KU.
    };

} // namespace LANGE

// TODO: provide an interface for different matrix types, i.e. GeMatrix,
//       TrMatrix, SbMatrix, HessenbergMatrix, vector types and scalars
//-- lascl ---------------------------------------------------------------------
template <typename IndexType, typename T, typename MA>
    typename RestrictTo<IsSame<typename MA::ElementType, T>::value, void>::Type
    lascl(LASCL::Type type, IndexType kl, IndexType ku,
          const T &cFrom, const T &cTo, MA &A);

template <typename IndexType, typename T, typename MA>
    typename RestrictTo<IsSame<MA, T>::value, void>::Type
    lascl(LASCL::Type type, IndexType kl, IndexType ku,
          const T &cFrom, const T &cTo, MA &A);

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename T, typename MA>
    typename RestrictTo<IsSame<typename MA::ElementType, T>::value, void>::Type
    lascl(LASCL::Type type, IndexType kl, IndexType ku,
          const T &cFrom, const T &cTo, MA &&A);

template <typename IndexType, typename T, typename VX>
    void
    lascl(LASCL::Type type, IndexType kl, IndexType ku,
          const T &cFrom, const T &cTo, DenseVector<VX> &x);

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_LASCL_H
