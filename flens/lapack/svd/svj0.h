/*
 *   Copyright (c) 2012, Michael Lehn
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
       SUBROUTINE DGSVJ0( JOBV, M, N, A, LDA, D, SVA, MV, V, LDV, EPS,
      $                   SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1)                                  --
 *
 *  -- Contributed by Zlatko Drmac of the University of Zagreb and     --
 *  -- Kresimir Veselic of the Fernuniversitaet Hagen                  --
 *  -- April 2011                                                      --
 *
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *
 * This routine is also part of SIGMA (version 1.23, October 23. 2008.)
 * SIGMA is a library of algorithms for highly accurate algorithms for
 * computation of SVD, PSVD, QSVD, (H,K)-SVD, and for solution of the
 * eigenvalue problems Hx = lambda M x, H M x = lambda x with H, M > 0.
 *
 */

#ifndef FLENS_LAPACK_SVD_SVJ0_H
#define FLENS_LAPACK_SVD_SVJ0_H 1

#include <flens/lapack/svd/svj.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== svj0 ======================================================================
template <typename MA, typename VD, typename VSVA, typename MV, typename VWORK>
    typename GeMatrix<MA>::IndexType
    svj0(SVJ::JobV                                  jobV,
         GeMatrix<MA>                               &A,
         DenseVector<VD>                            &d,
         DenseVector<VSVA>                          &sva,
         GeMatrix<MV>                               &V,
         const typename GeMatrix<MA>::ElementType   &eps,
         const typename GeMatrix<MA>::ElementType   &safeMin,
         const typename GeMatrix<MA>::ElementType   &tol,
         typename GeMatrix<MA>::IndexType           nSweep,
         DenseVector<VWORK>                         &work);

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VD, typename VSVA, typename MV, typename VWORK>
    typename MA::IndexType
    svj0(SVJ::JobV                                  jobV,
         MA                                         &&A,
         VD                                         &&d,
         VSVA                                       &&sva,
         MV                                         &&V,
         const typename MA::ElementType             &eps,
         const typename MA::ElementType             &safeMin,
         const typename MA::ElementType             &tol,
         typename MA::IndexType                     nSweep,
         VWORK                                      &&work);

} } // namespace lapack, flens

#endif // FLENS_LAPACK_SVD_SVJ0_H
