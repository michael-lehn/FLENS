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
       SUBROUTINE DGEJSV( JOBA, JOBU, JOBV, JOBR, JOBT, JOBP,
      $                   M, N, A, LDA, SVA, U, LDU, V, LDV,
      $                   WORK, LWORK, IWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1)                                    --
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

#ifndef FLENS_LAPACK_IMPL_JSV_H
#define FLENS_LAPACK_IMPL_JSV_H 1

#include <flens/lapack/typedefs.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

namespace JSV {

    // TODO: Find better (more expressive) names for these constants
    enum Accuracy  {
        accOptC = 'C',
        accOptE = 'E',
        accOptF = 'F',
        accOptG = 'G',
        accOptA = 'A',
        accOptR = 'R'
    };

    enum JobU {
        ComputeU   = 'U',
        FullsetU   = 'F',
        WorkspaceU = 'W',
        NoU        = 'N'
    };

    enum JobV {
        ComputeV        = 'V',
        JacobiRotationV = 'J',
        WorkspaceV      = 'W',
        NoV             = 'N'
    };
}

//== jsv =======================================================================
template <typename MA, typename VSVA, typename MU, typename MV,
          typename VWORK, typename VIWORK>
    typename GeMatrix<MA>::IndexType
    jsv(JSV::Accuracy             accuracy,
        JSV::JobU                 jobU,
        JSV::JobV                 jobV,
        bool                      restrictedRange,
        bool                      considerTransA,
        bool                      perturb,
        GeMatrix<MA>              &A,
        DenseVector<VSVA>         &sva,
        GeMatrix<MU>              &U,
        GeMatrix<MV>              &V,
        DenseVector<VWORK>        &work,
        DenseVector<VIWORK>       &iwork);

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VSVA, typename MU, typename MV,
          typename VWORK, typename VIWORK>
    typename MA::IndexType
    jsv(JSV::Accuracy             accuracy,
        JSV::JobU                 jobU,
        JSV::JobV                 jobV,
        bool                      restrictedRange,
        bool                      considerTransA,
        bool                      perturb,
        MA                        &&A,
        VSVA                      &&sva,
        MU                        &&U,
        MV                        &&V,
        VWORK                     &&work,
        VIWORK                    &&iwork);

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_JSV_H
