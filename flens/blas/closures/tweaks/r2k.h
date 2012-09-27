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

#ifndef FLENS_BLAS_CLOSURES_TWEAKS_R2K_H
#define FLENS_BLAS_CLOSURES_TWEAKS_R2K_H 1

#include <cxxblas/cxxblas.h>
#include <flens/blas/closures/tweaks/defaulteval.h>
#include <flens/blas/closures/tweaks/rk.h>
#include <flens/blas/operators/operators.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//
//  Real cases, i.e. not conjugated
//

//
//  Matrix closure of form:  a1*MA1*MB1^T + a2*MB2*MA2^T
//
template <typename MA1, typename MB1, typename MB2, typename MA2>
using MatrixClosureR2KU1_ =

    MatrixClosure<OpAdd,
                  MatrixClosureRKU1_<MA1, MB1>,
                  MatrixClosureRKU1_<MB2, MA2>
                 >;

//
//  Matrix closure of form:  MB + a1*MA1*MB1^T + a2*MB2*MA2^T
//
template <typename MB, typename MA1, typename MB1, typename MB2, typename MA2>
using MatrixClosureR2KU1 =

    MatrixClosure<OpAdd,
                  MatrixClosure<OpAdd,
                                MB,
                                MatrixClosureRKU1_<MA1, MB1>
                               >,
                  MatrixClosureRKU1_<MB2, MA2>
                 >;


//
//  Matrix closure of form:  MA1^T*MB1 + MB2^T*MA2
//
template <typename MA1, typename MB1, typename MB2, typename MA2>
using MatrixClosureR2KU2_ =

    MatrixClosure<OpAdd,
                  MatrixClosureRKU2_<MA1, MB1>,
                  MatrixClosureRKU2_<MB2, MA2>
                 >;

//
//  Matrix closure of form:  MB + MA1^T*MB1 + a2*MB2^T*MA2
//
template <typename MB, typename MA1, typename MB1, typename MB2, typename MA2>
using MatrixClosureR2KU2 =

    MatrixClosure<OpAdd,
                  MatrixClosure<OpAdd,
                                MB,
                                MatrixClosureRKU2_<MA1, MB1>
                               >,
                  MatrixClosureRKU2_<MB2, MA2>
                 >;


//
//  Matrix closure of form:  a1*MA1^T*MB1 + a2*MB2^T*MA2
//
template <typename SV, typename MA1, typename MB1, typename MB2, typename MA2>
using MatrixClosureR2KU3_ =

    MatrixClosure<OpAdd,
                  MatrixClosureRKU3_<SV, MA1, MB1>,
                  MatrixClosureRKU3_<SV, MB2, MA2>
                 >;

//
//  Matrix closure of form:  MB + MA1^T*MB1 + a2*MB2^T*MA2
//
template <typename MB, typename SV,
          typename MA1, typename MB1, typename MB2, typename MA2>
using MatrixClosureR2KU3 =

    MatrixClosure<OpAdd,
                  MatrixClosure<OpAdd,
                                MB,
                                MatrixClosureRKU3_<SV, MA1, MB1>
                               >,
                  MatrixClosureRKU3_<SV, MB2, MA2>
                 >;

//
//  Complex cases
//

//
//  Matrix closure of form:  a1*MA1*MB1^H + a2*MB2*MA2^H
//
template <typename MA1, typename MB1, typename MB2, typename MA2>
using MatrixClosureR2KC1_ =

    MatrixClosure<OpAdd,
                  MatrixClosureRKC1_<MA1, MB1>,
                  MatrixClosureRKC1_<MB2, MA2>
                 >;

//
//  Matrix closure of form:  MB + a1*MA1*MB1^H + a2*MB2*MA2^H
//
template <typename MB, typename MA1, typename MB1, typename MB2, typename MA2>
using MatrixClosureR2KC1 =

    MatrixClosure<OpAdd,
                  MatrixClosure<OpAdd,
                                MB,
                                MatrixClosureRKC1_<MA1, MB1>
                               >,
                  MatrixClosureRKC1_<MB2, MA2>
                 >;


//
//  Matrix closure of form:  MA1^H*MB1 + MB2^H*MA2
//
template <typename MA1, typename MB1, typename MB2, typename MA2>
using MatrixClosureR2KC2_ =

    MatrixClosure<OpAdd,
                  MatrixClosureRKC2_<MA1, MB1>,
                  MatrixClosureRKC2_<MB2, MA2>
                 >;

//
//  Matrix closure of form:  MB + MA1^H*MB1 + a2*MB2^H*MA2
//
template <typename MB, typename MA1, typename MB1, typename MB2, typename MA2>
using MatrixClosureR2KC2 =

    MatrixClosure<OpAdd,
                  MatrixClosure<OpAdd,
                                MB,
                                MatrixClosureRKC2_<MA1, MB1>
                               >,
                  MatrixClosureRKC2_<MB2, MA2>
                 >;


//
//  Matrix closure of form:  a1*MA1^H*MB1 + a2*MB2^H*MA2
//
template <typename SV, typename MA1, typename MB1, typename MB2, typename MA2>
using MatrixClosureR2KC3_ =

    MatrixClosure<OpAdd,
                  MatrixClosureRKC3_<SV, MA1, MB1>,
                  MatrixClosureRKC3_<SV, MB2, MA2>
                 >;

//
//  Matrix closure of form:  MB + MA1^H*MB1 + a2*MB2^H*MA2
//
template <typename MB, typename SV,
          typename MA1, typename MB1, typename MB2, typename MA2>
using MatrixClosureR2KC3 =

    MatrixClosure<OpAdd,
                  MatrixClosure<OpAdd,
                                MB,
                                MatrixClosureRKC3_<SV, MA1, MB1>
                               >,
                  MatrixClosureRKC3_<SV, MB2, MA2>
                 >;


//-- SymmetricMatrix -----------------------------------------------------------

//
//  Case 1: trans==NoTrans
//

//
//  C += a1*A1*B1^T + a2*B2*A2^T
//
template <typename ALPHA,
          typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KU1_<MA1, MB1,
                                                        MB2, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureR2KU1_<MA1, MB1, MB2, MA2> &aABt_aBAt,
         MC &C);

//
//  C = a*A*B^T + a*B*A^T
//
template <typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KU1_<MA1, MB1,
                                                        MB2, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureR2KU1_<MA1, MB1, MB2, MA2> &aABt_aBAt,
         MC &C);

//
//  C = b*D + a*A*B^T + a*B*A^T
//
template <typename MD,
          typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KU1<MD, MA1, MB1,
                                                           MB2, MA2> >::value
                     && IsMatrix<MD>::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureR2KU1<MD, MA1, MB1, MB2, MA2> &bD_aABt_aBAt,
         MC &C);

//
//  Case 2: trans==Trans, alpha==1
//

//
//  C += A^T*B + B^T*A
//
template <typename ALPHA,
          typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KU2_<MA1, MB1,
                                                        MB2, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureR2KU2_<MA1, MB1, MB2, MA2> &AtB_BtA,
         MC &C);

//
//  C = A^T*B + B^T*A
//
template <typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KU2_<MA1, MB1,
                                                        MB2, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureR2KU2_<MA1, MB1, MB2, MA2> &AtB_BtA,
         MC &C);

//
//  C = b*D + A^T*B + B^T*A
//
template <typename MD,
          typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KU2<MD, MA1, MB1,
                                                           MB2, MA2> >::value
                     && IsMatrix<MD>::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureR2KU2<MD, MA1, MB1, MB2, MA2> &bD_AtB_BtA,
         MC &C);

//
//  Case 3: trans==Trans, alpha
//

//
//  C += a*A^T*B + a*B^T*A
//
template <typename ALPHA,
          typename SV, typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KU3_<SV, MA1, MB1,
                                                            MB2, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureR2KU3_<SV, MA1, MB1, MB2, MA2> &aAtB_aBtA,
         MC &C);

//
//  C = a*A^T*B + a*B^T*A
//
template <typename SV, typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KU3_<SV, MA1, MB1,
                                                            MB2, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureR2KU3_<SV, MA1, MB1, MB2, MA2> &aAtB_aBtA,
         MC &C);

//
//  C = b*D + a*A^T*B + a*B^T*A
//
template <typename MD,
          typename SV, typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KU3<MD, SV,
                                                       MA1, MB1,
                                                       MB2, MA2> >::value
                     && IsMatrix<MD>::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureR2KU3<MD, SV, MA1, MB1, MB2, MA2> &bD_aAtB_aBtA,
         MC &C);


//-- HermitianMatrix -----------------------------------------------------------

//
//  Case 1: trans==NoTrans
//

//
//  C += a1*A1*B1^H + a2*B2*A2^H
//
template <typename ALPHA,
          typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KC1_<MA1, MB1,
                                                        MB2, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureR2KC1_<MA1, MB1, MB2, MA2> &aABt_aBAt,
         MC &C);

//
//  C = a*A*B^H + a*B*A^H
//
template <typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KC1_<MA1, MB1,
                                                        MB2, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureR2KC1_<MA1, MB1, MB2, MA2> &aABt_aBAt,
         MC &C);

//
//  C = b*D + a*A*B^H + a*B*A^H
//
template <typename MD,
          typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KC1<MD, MA1, MB1,
                                                           MB2, MA2> >::value
                     && IsMatrix<MD>::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureR2KC1<MD, MA1, MB1, MB2, MA2> &bD_aABt_aBAt,
         MC &C);

//
//  Case 2: trans==ConjTrans, alpha==1
//

//
//  C += A^H*B + B^H*A
//
template <typename ALPHA,
          typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KC2_<MA1, MB1,
                                                        MB2, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureR2KC2_<MA1, MB1, MB2, MA2> &AtB_BtA,
         MC &C);

//
//  C = A^H*B + B^H*A
//
template <typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KC2_<MA1, MB1,
                                                        MB2, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureR2KC2_<MA1, MB1, MB2, MA2> &AtB_BtA,
         MC &C);

//
//  C = b*D + A^H*B + B^H*A
//
template <typename MD,
          typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KC2<MD, MA1, MB1,
                                                           MB2, MA2> >::value
                     && IsMatrix<MD>::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureR2KC2<MD, MA1, MB1, MB2, MA2> &bD_AtB_BtA,
         MC &C);

//
//  Case 3: trans==ConjTrans, alpha
//

//
//  C += a*A^H*B + a*B^H*A
//
template <typename ALPHA,
          typename SV, typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KC3_<SV, MA1, MB1,
                                                            MB2, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureR2KC3_<SV, MA1, MB1, MB2, MA2> &aAtB_aBtA,
         MC &C);

//
//  C = a*A^H*B + a*B^H*A
//
template <typename SV, typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KC3_<SV, MA1, MB1,
                                                            MB2, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureR2KC3_<SV, MA1, MB1, MB2, MA2> &aAtB_aBtA,
         MC &C);

//
//  C = b*D + a*A^H*B + a*B^H*A
//
template <typename MD,
          typename SV, typename MA1, typename MB1, typename MB2, typename MA2,
          typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureR2KC3<MD, SV,
                                                       MA1, MB1,
                                                       MB2, MA2> >::value
                     && IsMatrix<MD>::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MB1>::value
                     && IsMatrix<MB2>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureR2KC3<MD, SV, MA1, MB1, MB2, MA2> &bD_aAtB_aBtA,
         MC &C);


} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_TWEAKS_R2K_H
