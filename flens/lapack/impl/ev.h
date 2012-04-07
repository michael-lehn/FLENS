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
      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
     $                  LDVR, WORK, LWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_EV_H
#define FLENS_LAPACK_IMPL_EV_H 1

#include <flens/aux/aux.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== ev_wsq (workspace query) ==================================================
template <typename MA>
    Pair<typename GeMatrix<MA>::IndexType>
    ev_wsq(bool computeVL, bool computeVR, GeMatrix<MA> &A);

//== ev ========================================================================
template <typename MA, typename VWR, typename VWI, typename MVL, typename MVR,
          typename VWORK>
    typename GeMatrix<MA>::IndexType
    ev(bool computeVL, bool computeVR,
       GeMatrix<MA> &A,
       DenseVector<VWR> &wr, DenseVector<VWI> &wi,
       GeMatrix<MVL> &VL, GeMatrix<MVR> &VR,
       DenseVector<VWORK> &work);

//-- forwarding ----------------------------------------------------------------
template <typename MA>
    typename MA::IndexType
    ev_wsq(bool computeVL, bool computeVR, MA &&A);

template <typename MA, typename VWR, typename VWI, typename MVL, typename MVR,
          typename VWORK>
    typename MA::IndexType
    ev(bool computeVL, bool computeVR,
       MA &&A,  VWR &&wr, VWI &&wi,
       MVL &&VL, MVR &&VR,
       VWORK &&work);

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_EV_H
