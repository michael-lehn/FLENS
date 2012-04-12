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
      SUBROUTINE DTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
     $                   LDVR, S, SEP, MM, M, WORK, LDWORK, IWORK,
     $                   INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_TRSNA_H
#define FLENS_LAPACK_IMPL_TRSNA_H 1

#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

namespace TRSNA {

    enum Job {
        EigenvaluesOnly  = 'E', // = 'E': for eigenvalues only (S);
        EigenvectorsOnly = 'V', // = 'V': for eigenvectors only (SEP);
        Both             = 'B'  // = 'B': for both eigenvalues and
                                //        eigenvectors (S and SEP).
    };

    enum HowMany {
        All      = 'A', // = 'A': compute condition numbers for all eigenpairs
        Selected = 'S'  // = 'S': compute condition numbers for selected
                        //        eigenpairs specified by the array SELECT.
    };
}

//== trsna =====================================================================
template <typename VSELECT, typename MT, typename MVL, typename MVR,
          typename VS, typename VSEP, typename MM, typename M,
          typename MWORK, typename VIWORK>
    void
    trsna(TRSNA::Job                    job,
          TRSNA::HowMany                howMany,
          const DenseVector<VSELECT>    &select,
          const GeMatrix<MT>            &T,
          const GeMatrix<MVL>           &VL,
          const GeMatrix<MVR>           &VR,
          DenseVector<VS>               &s,
          DenseVector<VSEP>             &sep,
          const MM                      &mm,
          M                             &m,
          GeMatrix<MWORK>               &work,
          DenseVector<VIWORK>           &iWork);

//-- forwarding ----------------------------------------------------------------
template <typename VSELECT, typename MT, typename MVL, typename MVR,
          typename VS, typename VSEP, typename MM, typename M,
          typename MWORK, typename VIWORK>
    void
    trsna(TRSNA::Job                    job,
          TRSNA::HowMany                howMany,
          const VSELECT                 &select,
          const MT                      &T,
          const MVL                     &VL,
          const MVR                     &VR,
          VS                            &&s,
          VSEP                          &&sep,
          const MM                      &mm,
          M                             &&m,
          MWORK                         &&work,
          VIWORK                        &&iWork);

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_TRSNA_H
