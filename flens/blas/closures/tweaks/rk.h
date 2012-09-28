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

#ifndef FLENS_BLAS_CLOSURES_TWEAKS_RK_H
#define FLENS_BLAS_CLOSURES_TWEAKS_RK_H 1

#include <cxxblas/cxxblas.h>
#include <flens/blas/closures/tweaks/defaulteval.h>
#include <flens/blas/closures/tweaks/r.h>
#include <flens/blas/operators/operators.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//
//  Real cases, i.e. not conjugated
//

//
//  Matrix closure of form:  MA1*MA2^T
//
template <typename MA1, typename MA2>
using MatrixClosureRKU1_ =

    MatrixClosure<OpMult,
                  MA1,
                  MatrixClosureOpTrans<MA2>
                 >;

//
//  Matrix closure of form:  MB + MA1*MA2^T
//
template <typename MB, typename MA1, typename MA2>
using MatrixClosureRKU1 =

    MatrixClosure<OpAdd,
                  MB,
                  MatrixClosureRKU1_<MA1, MA2>
                 >;

//
//  Matrix closure of form:  MA1^T*MA2
//
template <typename MA1, typename MA2>
using MatrixClosureRKU2_ =

    MatrixClosure<OpMult,
                  MatrixClosureOpTrans<MA1>,
                  MA2
                 >;

//
//  Matrix closure of form:  MB + MA1^T*MA2
//
template <typename MB, typename MA1, typename MA2>
using MatrixClosureRKU2 =

    MatrixClosure<OpAdd,
                  MB,
                  MatrixClosureRKU2_<MA1, MA2>
                 >;

//
//  Matrix closure of form:  (SV*MA1^T)*MA2
//
template <typename SV, typename MA1, typename MA2>
using MatrixClosureRKU3_ =

    MatrixClosure<OpMult,
                  MatrixClosure<OpMult,
                                ScalarValue<SV>,
                                MatrixClosureOpTrans<MA1>
                               >,
                  MA2
                 >;

//
//  Matrix closure of form:  MB + (SV*MA1^T)*MA2
//
template <typename MB, typename SV, typename MA1, typename MA2>
using MatrixClosureRKU3 =

    MatrixClosure<OpAdd,
                  MB,
                  MatrixClosureRKU3_<SV, MA1, MA2>
                 >;


//
//  Complex cases
//

//
//  Matrix closure of form:  MA1*MA2^H
//
template <typename MA1, typename MA2>
using MatrixClosureRKC1_ =

    MatrixClosure<OpMult,
                  MA1,
                  MatrixClosureOpConjTrans<MA2>
                 >;

//
//  Matrix closure of form:  MB + MA1*MA2^H
//
template <typename MB, typename MA1, typename MA2>
using MatrixClosureRKC1 =

    MatrixClosure<OpAdd,
                  MB,
                  MatrixClosureRKC1_<MA1, MA2>
                 >;

//
//  Matrix closure of form:  MA1^H*MA2
//
template <typename MA1, typename MA2>
using MatrixClosureRKC2_ =

    MatrixClosure<OpMult,
                  MatrixClosureOpConjTrans<MA1>,
                  MA2
                 >;

//
//  Matrix closure of form:  MB + MA1^H*MA2
//
template <typename MB, typename MA1, typename MA2>
using MatrixClosureRKC2 =

    MatrixClosure<OpAdd,
                  MB,
                  MatrixClosureRKC2_<MA1, MA2>
                 >;

//
//  Matrix closure of form:  (SV*MA1^H)*MA2
//
template <typename SV, typename MA1, typename MA2>
using MatrixClosureRKC3_ =

    MatrixClosure<OpMult,
                  MatrixClosure<OpMult,
                                ScalarValue<SV>,
                                MatrixClosureOpConjTrans<MA1>
                               >,
                  MA2
                 >;

//
//  Matrix closure of form:  MB + (SV*MA1^H)*MA2
//
template <typename MB, typename SV, typename MA1, typename MA2>
using MatrixClosureRKC3 =

    MatrixClosure<OpAdd,
                  MB,
                  MatrixClosureRKC3_<SV, MA1, MA2>
                 >;



//-- SymmetricMatrix -----------------------------------------------------------

//
//  Case 1: trans==NoTrans
//

//
//  C += a*A*A^T
//
template <typename ALPHA, typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKU1_<MA1, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureRKU1_<MA1, MA2> &aAAt,
         MC &C);

//
//  C = a*A*A^T
//
template <typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKU1_<MA1, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureRKU1_<MA1, MA2> &aAAt,
         MC &C);

//
//  C = b*B + a*A*A^T
//
template <typename MB, typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKU1<MB, MA1, MA2> >::value
                     && IsMatrix<MB>::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureRKU1<MB, MA1, MA2> &bB_aAAt,
         MC &C);

//
//  Case 2: trans==Trans, alpha==1
//

//
//  C += A^T*A
//
template <typename ALPHA, typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKU2_<MA1, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureRKU2_<MA1, MA2> &AtA,
         MC &C);

//
//  C = A^T*A
//
template <typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKU2_<MA1, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureRKU2_<MA1, MA2> &AtA,
         MC &C);

//
//  C = b*B + A^T*A
//
template <typename MB, typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKU2<MB, MA1, MA2> >::value
                     && IsMatrix<MB>::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureRKU2<MB, MA1, MA2> &bB_AtA,
         MC &C);

//
//  Case 3: trans==Trans, alpha
//

//
//  C += a*A^T*A
//
template <typename ALPHA, typename SV, typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKU3_<SV, MA1, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureRKU3_<SV, MA1, MA2> &aAtA,
         MC &C);

//
//  C = a*A^T*A
//
template <typename SV, typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKU3_<SV, MA1, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureRKU3_<SV, MA1, MA2> &aAtA,
         MC &C);

//
//  C = b*B + a*A^T*A
//
template <typename MB, typename SV, typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKU3<MB, SV, MA1, MA2> >::value
                     && IsMatrix<MB>::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsSymmetricMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureRKU3<MB, SV, MA1, MA2> &bB_AtA,
         MC &C);


//-- HermitianMatrix -----------------------------------------------------------

//
//  Case 1: trans==NoTrans
//

//
//  C += a*A*A^H
//
template <typename ALPHA, typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKC1_<MA1, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureRKC1_<MA1, MA2> &aAAh,
         MC &C);

//
//  C = a*A*A^H
//
template <typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKC1_<MA1, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureRKC1_<MA1, MA2> &aAAh,
         MC &C);

//
//  C = b*B + a*A*A^H
//
template <typename MB, typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKC1<MB, MA1, MA2> >::value
                     && IsMatrix<MB>::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureRKC1<MB, MA1, MA2> &bB_aAAh,
         MC &C);

//
//  Case 2: trans==Trans, alpha==1
//

//
//  C += A^H*A
//
template <typename ALPHA, typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKC2_<MA1, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureRKC2_<MA1, MA2> &AhA,
         MC &C);

//
//  C = A^H*A
//
template <typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKC2_<MA1, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureRKC2_<MA1, MA2> &AhA,
         MC &C);


//
//  C = b*B + A^H*A
//
template <typename MB, typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKC2<MB, MA1, MA2> >::value
                     && IsMatrix<MB>::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureRKC2<MB, MA1, MA2> &bB_AhA,
         MC &C);

//
//  Case 3: trans==Trans, alpha
//

//
//  C += a*A^H*A
//
template <typename ALPHA, typename SV, typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKC3_<SV, MA1, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureRKC3_<SV, MA1, MA2> &aAtA,
         MC &C);

//
//  C = a*A^H*A
//
template <typename SV, typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKC3_<SV, MA1, MA2> >::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureRKC3_<SV, MA1, MA2> &aAtA,
         MC &C);

//
//  C = b*B + a*A^H*A
//
template <typename MB, typename SV, typename MA1, typename MA2, typename MC>
    typename RestrictTo<DefaultEval<MatrixClosureRKC3<MB, SV, MA1, MA2> >::value
                     && IsMatrix<MB>::value
                     && IsMatrix<MA1>::value
                     && IsMatrix<MA2>::value
                     && IsHermitianMatrix<MC>::value,
             void>::Type
    copy(Transpose trans,
         const MatrixClosureRKC3<MB, SV, MA1, MA2> &bB_AtA,
         MC &C);


} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_TWEAKS_RK_H
