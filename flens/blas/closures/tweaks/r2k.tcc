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

#ifndef FLENS_BLAS_CLOSURES_TWEAKS_R2K_TCC
#define FLENS_BLAS_CLOSURES_TWEAKS_R2K_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/closures/closures.h>
#include <flens/blas/level1/level1.h>
#include <flens/blas/level2/level2.h>
#include <flens/blas/level3/level3.h>
#include <flens/typedefs.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//-- SymmetricMatrix -----------------------------------------------------------

//
//  Case1: trans==NoTrans
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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MA1>::Remainder    RMA;
    typedef typename PruneScaling<MA1>::ScalingType  SMA;

    typedef typename PruneScaling<MB2>::Remainder    RMB;
    typedef typename PruneScaling<MB2>::ScalingType  SMB;

    const SMA &a1  = PruneScaling<MA1>::getFactor(aABt_aBAt.left().left());
    const RMA &A1  = PruneScaling<MA1>::getRemainder(aABt_aBAt.left().left());
    const auto &B1 = aABt_aBAt.left().right().right();

    const SMB &a2  = PruneScaling<MB2>::getFactor(aABt_aBAt.right().left());
    const RMB &B2  = PruneScaling<MB2>::getRemainder(aABt_aBAt.right().left());
    const auto &A2 = aABt_aBAt.right().right().right();

    if (identical(A1, A2) && identical(B1, B2) && (a1==a2)) {
        blas::r2k(trans, alpha*a1, A1, B1, SMA(1), C.impl());
        return;
    } else {
        blas::axpy(trans, alpha, aABt_aBAt.left(), C.impl());
        blas::axpy(trans, alpha, aABt_aBAt.right(), C.impl());
    }
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MA1>::Remainder    RMA;
    typedef typename PruneScaling<MA1>::ScalingType  SMA;

    typedef typename PruneScaling<MB2>::Remainder    RMB;
    typedef typename PruneScaling<MB2>::ScalingType  SMB;

    const SMA &a1  = PruneScaling<MA1>::getFactor(aABt_aBAt.left().left());
    const RMA &A1  = PruneScaling<MA1>::getRemainder(aABt_aBAt.left().left());
    const auto &B1 = aABt_aBAt.left().right().right();

    const SMB &a2  = PruneScaling<MB2>::getFactor(aABt_aBAt.right().left());
    const RMB &B2  = PruneScaling<MB2>::getRemainder(aABt_aBAt.right().left());
    const auto &A2 = aABt_aBAt.right().right().right();

    if (identical(A1, A2) && identical(B1, B2) && (a1==a2)) {
        blas::r2k(trans, a1, A1, B1, SMA(0), C.impl());
        return;
    } else {
        blas::copy(trans, aABt_aBAt.left(), C.impl());
        blas::axpy(trans, SMB(1), aABt_aBAt.right(), C.impl());
    }
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MD>::Remainder     RMD;
    typedef typename PruneScaling<MD>::ScalingType   SMD;

    typedef typename PruneScaling<MA1>::Remainder    RMA;
    typedef typename PruneScaling<MA1>::ScalingType  SMA;

    typedef typename PruneScaling<MB2>::Remainder    RMB;
    typedef typename PruneScaling<MB2>::ScalingType  SMB;

    const auto &bD   = bD_aABt_aBAt.left().left();
    const auto &aABt = bD_aABt_aBAt.left().right();
    const auto &aBAt = bD_aABt_aBAt.right();

    const SMD &b   = PruneScaling<MD>::getFactor(bD);
    const RMD &D   = PruneScaling<MD>::getRemainder(bD);

    const SMA &a1  = PruneScaling<MA1>::getFactor(aABt.left());
    const RMA &A1  = PruneScaling<MA1>::getRemainder(aABt.left());
    const auto &B1 = aABt.right().right();

    const SMB &a2  = PruneScaling<MB2>::getFactor(aBAt.left());
    const RMB &B2  = PruneScaling<MB2>::getRemainder(aBAt.left());
    const auto &A2 = aBAt.right().right();

    if (identical(D, C.impl())
     && identical(A1, A2)
     && identical(B1, B2)
     && (a1==a2))
    {
        blas::r2k(trans, a1, A1, B1, b, C.impl());
        return;
    } else {
        blas::copy(trans, D, C.impl());
        blas::scal(b, C.impl());
        blas::axpy(trans, SMB(1), aABt, C.impl());
        blas::axpy(trans, SMB(1), aBAt, C.impl());
    }
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    const auto &A1 = AtB_BtA.left().left().left();
    const auto &B1 = AtB_BtA.left().right();

    const auto &B2 = AtB_BtA.right().left().left();
    const auto &A2 = AtB_BtA.right().right();

    if (identical(A1, A2) && identical(B1, B2)) {
        trans = Transpose(trans^Trans);
        blas::r2k(trans, alpha, A1, B1, ALPHA(1), C.impl());
        return;
    } else {
        std::cerr << "warning" << std::endl;
        blas::axpy(trans, alpha, AtB_BtA.left(), C.impl());
        blas::axpy(trans, alpha, AtB_BtA.right(), C.impl());
    }
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    const auto &A1 = AtB_BtA.left().left().left();
    const auto &B1 = AtB_BtA.left().right();

    const auto &B2 = AtB_BtA.right().left().left();
    const auto &A2 = AtB_BtA.right().right();

    typedef typename MC::Impl::ElementType  ALPHA;

    if (identical(A1, A2) && identical(B1, B2)) {
        trans = Transpose(trans^Trans);
        blas::r2k(trans, ALPHA(1), A1, B1, ALPHA(0), C.impl());
        return;
    } else {

        std::cerr << "warning" << std::endl;
        blas::copy(trans, AtB_BtA.left(), C.impl());
        blas::axpy(trans, ALPHA(1), AtB_BtA.right(), C.impl());
    }
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MD>::Remainder     RMD;
    typedef typename PruneScaling<MD>::ScalingType   SMD;

    const auto &bD  = bD_AtB_BtA.left().left();
    const auto &AtB = bD_AtB_BtA.left().right();
    const auto &BtA = bD_AtB_BtA.right();

    const SMD &b    = PruneScaling<MD>::getFactor(bD);
    const RMD &D    = PruneScaling<MD>::getRemainder(bD);

    const auto &A1  = AtB.left().left();
    const auto &B1  = AtB.right();

    const auto &B2  = BtA.left().left();
    const auto &A2  = BtA.right();

    typedef typename MC::Impl::ElementType  ALPHA;

    if (identical(D, C.impl()) && identical(A1, A2) && identical(B1, B2)) {
        trans = Transpose(trans^Trans);
        blas::r2k(trans, ALPHA(1), A1, B1, b, C.impl());
        return;
    } else {
        blas::copy(trans, D, C.impl());
        blas::scal(b, C.impl());
        blas::axpy(trans, SMD(1), AtB, C.impl());
        blas::axpy(trans, SMD(1), BtA, C.impl());
    }
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    const auto &aAtB = aAtB_aBtA.left();
    const auto &aBtA = aAtB_aBtA.right();

    const auto &a1 = aAtB.left().left().value();
    const auto &A1 = aAtB.left().right().right();
    const auto &B1 = aAtB.right();

    const auto &a2 = aBtA.left().left().value();
    const auto &B2 = aBtA.left().right().right();
    const auto &A2 = aBtA.right();

    if ((a1==a2) && identical(A1, A2) && identical(B1, B2)) {
        trans = Transpose(trans^Trans);
        blas::r2k(trans, alpha*a1, A1, B1, ALPHA(1), C.impl());
        return;
    } else {
        std::cerr << "warning" << std::endl;
        blas::axpy(trans, alpha, aAtB, C.impl());
        blas::axpy(trans, alpha, aBtA, C.impl());
    }
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    const auto &aAtB = aAtB_aBtA.left();
    const auto &aBtA = aAtB_aBtA.right();

    const auto &a1 = aAtB.left().left().value();
    const auto &A1 = aAtB.left().right().right();
    const auto &B1 = aAtB.right();

    const auto &a2 = aBtA.left().left().value();
    const auto &B2 = aBtA.left().right().right();
    const auto &A2 = aBtA.right();

    typedef typename MC::Impl::ElementType  ALPHA;

    if ((a1==a2) && identical(A1, A2) && identical(B1, B2)) {
        trans = Transpose(trans^Trans);
        blas::r2k(trans, a1, A1, B1, ALPHA(0), C.impl());
        return;
    } else {
        std::cerr << "warning" << std::endl;
        blas::copy(trans, aAtB, C.impl());
        blas::axpy(trans, ALPHA(1), aBtA, C.impl());
    }
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MD>::Remainder     RMD;
    typedef typename PruneScaling<MD>::ScalingType   SMD;

    const auto &bD   = bD_aAtB_aBtA.left().left();
    const auto &aAtB = bD_aAtB_aBtA.left().right();
    const auto &aBtA = bD_aAtB_aBtA.right();

    const SMD &b    = PruneScaling<MD>::getFactor(bD);
    const RMD &D    = PruneScaling<MD>::getRemainder(bD);

    const auto &a1 = aAtB.left().left().value();
    const auto &A1 = aAtB.left().right().right();
    const auto &B1 = aAtB.right();

    const auto &a2 = aBtA.left().left().value();
    const auto &B2 = aBtA.left().right().right();
    const auto &A2 = aBtA.right();

    typedef typename MC::Impl::ElementType  ALPHA;

    if ((a1==a2) && identical(A1, A2) && identical(B1, B2)) {
        trans = Transpose(trans^Trans);
        blas::r2k(trans, a1, A1, B1, b, C.impl());
        return;
    } else {
        blas::copy(trans, D, C.impl());
        blas::scal(b, C.impl());
        blas::axpy(trans, ALPHA(1), aAtB, C.impl());
        blas::axpy(trans, ALPHA(1), aBtA, C.impl());
    }
}

//-- HermitianMatrix -----------------------------------------------------------

//
//  Case1: trans==NoTrans
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
     const MatrixClosureR2KC1_<MA1, MB1, MB2, MA2> &aABh_aBAh,
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MA1>::Remainder    RMA;
    typedef typename PruneScaling<MA1>::ScalingType  SMA;

    typedef typename PruneScaling<MB2>::Remainder    RMB;
    typedef typename PruneScaling<MB2>::ScalingType  SMB;

    const auto &aABh = aABh_aBAh.left();
    const auto &aBAh = aABh_aBAh.right();

    const SMA &a1  = PruneScaling<MA1>::getFactor(aABh.left());
    const RMA &A1  = PruneScaling<MA1>::getRemainder(aABh.left());
    const auto &B1 = aABh.right().right().right();

    const SMB &a2  = PruneScaling<MB2>::getFactor(aBAh.left());
    const RMB &B2  = PruneScaling<MB2>::getRemainder(aBAh.left());
    const auto &A2 = aBAh.right().right().right();

    if (identical(A1, A2)
     && identical(B1, B2)
     && (a1==cxxblas::conjugate(a2)))
    {
        blas::r2k(trans, alpha*a1, A1, B1, SMA(1), C.impl());
        return;
    } else {
        std::cerr << "warning" << std::endl;
        blas::axpy(trans, alpha, aABh, C.impl());
        blas::axpy(trans, alpha, aBAh, C.impl());
    }
}

//
//  C = a1*A1*B1^H + a2*B2*A2^H
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
     const MatrixClosureR2KC1_<MA1, MB1, MB2, MA2> &aABh_aBAh,
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MA1>::Remainder    RMA;
    typedef typename PruneScaling<MA1>::ScalingType  SMA;

    typedef typename PruneScaling<MB2>::Remainder    RMB;
    typedef typename PruneScaling<MB2>::ScalingType  SMB;

    const auto &aABh = aABh_aBAh.left();
    const auto &aBAh = aABh_aBAh.right();

    const SMA &a1  = PruneScaling<MA1>::getFactor(aABh.left());
    const RMA &A1  = PruneScaling<MA1>::getRemainder(aABh.left());
    const auto &B1 = aABh.right().right().right();

    const SMB &a2  = PruneScaling<MB2>::getFactor(aBAh.left());
    const RMB &B2  = PruneScaling<MB2>::getRemainder(aBAh.left());
    const auto &A2 = aBAh.right().right().right();

    if (identical(A1, A2)
     && identical(B1, B2)
     && (a1==cxxblas::conjugate(a2)))
    {
        blas::r2k(trans, a1, A1, B1, SMA(0), C.impl());
        return;
    } else {
        blas::copy(trans, aABh, C.impl());
        blas::axpy(trans, SMB(1), aBAh, C.impl());
    }
}

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
     const MatrixClosureR2KC1<MD, MA1, MB1, MB2, MA2> &bD_aABh_aBAh,
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MD>::Remainder     RMD;
    typedef typename PruneScaling<MD>::ScalingType   SMD;

    typedef typename PruneScaling<MA1>::Remainder    RMA;
    typedef typename PruneScaling<MA1>::ScalingType  SMA;

    typedef typename PruneScaling<MB2>::Remainder    RMB;
    typedef typename PruneScaling<MB2>::ScalingType  SMB;

    const auto &bD   = bD_aABh_aBAh.left().left();
    const auto &aABh = bD_aABh_aBAh.left().right();
    const auto &aBAh = bD_aABh_aBAh.right();

    const SMD &b   = PruneScaling<MD>::getFactor(bD);
    const RMD &D   = PruneScaling<MD>::getRemainder(bD);

    const SMA &a1  = PruneScaling<MA1>::getFactor(aABh.left());
    const RMA &A1  = PruneScaling<MA1>::getRemainder(aABh.left());
    const auto &B1 = aABh.right().right().right();

    const SMB &a2  = PruneScaling<MB2>::getFactor(aBAh.left());
    const RMB &B2  = PruneScaling<MB2>::getRemainder(aBAh.left());
    const auto &A2 = aBAh.right().right().right();

    if (identical(D, C.impl())
     && identical(A1, A2)
     && identical(B1, B2)
     && (a1==cxxblas::conjugate(a2)))
    {
        blas::r2k(trans, a1, A1, B1, b, C.impl());
        return;
    } else {
        blas::copy(trans, D, C.impl());
        blas::scal(b, C.impl());
        blas::axpy(trans, SMB(1), aABh, C.impl());
        blas::axpy(trans, SMB(1), aBAh, C.impl());
    }
}

//
//  Case 2: trans==Trans, alpha==1
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
     const MatrixClosureR2KC2_<MA1, MB1, MB2, MA2> &AhB_BhA,
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    const auto &AhB = AhB_BhA.left();
    const auto &BhA = AhB_BhA.right();

    const auto &A1 = AhB.left().left().left();
    const auto &B1 = AhB.right();

    const auto &B2 = BhA.left().left().left();
    const auto &A2 = BhA.right();

    if (identical(A1, A2) && identical(B1, B2)) {
        trans = Transpose(trans^ConjTrans);
        blas::r2k(trans, alpha, A1, B1, ALPHA(1), C.impl());
        return;
    } else {
        std::cerr << "warning" << std::endl;
        blas::axpy(trans, alpha, AhB, C.impl());
        blas::axpy(trans, alpha, BhA, C.impl());
    }
}

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
     const MatrixClosureR2KC2_<MA1, MB1, MB2, MA2> &AhB_BhA,
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    const auto &AhB = AhB_BhA.left();
    const auto &BhA = AhB_BhA.right();

    const auto &A1 = AhB.left().left().left();
    const auto &B1 = AhB.right();

    const auto &B2 = BhA.left().left().left();
    const auto &A2 = BhA.right();

    typedef typename MC::Impl::ElementType  ALPHA;

    if (identical(A1, A2) && identical(B1, B2)) {
        trans = Transpose(trans^ConjTrans);
        blas::r2k(trans, ALPHA(1), A1, B1, ALPHA(0), C.impl());
        return;
    } else {

        std::cerr << "warning" << std::endl;
        blas::copy(trans, AhB, C.impl());
        blas::axpy(trans, ALPHA(1), BhA, C.impl());
    }
}

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
     const MatrixClosureR2KC2<MD, MA1, MB1, MB2, MA2> &bD_AhB_BhA,
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MD>::Remainder     RMD;
    typedef typename PruneScaling<MD>::ScalingType   SMD;

    const auto &bD  = bD_AhB_BhA.left().left();
    const auto &AhB = bD_AhB_BhA.left().right();
    const auto &BhA = bD_AhB_BhA.right();

    const SMD &b    = PruneScaling<MD>::getFactor(bD);
    const RMD &D    = PruneScaling<MD>::getRemainder(bD);

    const auto &A1  = AhB.left().left().left();
    const auto &B1  = AhB.right();

    const auto &B2  = BhA.left().left().left();
    const auto &A2  = BhA.right();

    typedef typename MC::Impl::ElementType  ALPHA;

    if (identical(D, C.impl()) && identical(A1, A2) && identical(B1, B2)) {
        trans = Transpose(trans^ConjTrans);
        blas::r2k(trans, ALPHA(1), A1, B1, b, C.impl());
        return;
    } else {
        blas::copy(trans, D, C.impl());
        blas::scal(b, C.impl());
        blas::axpy(trans, SMD(1), AhB, C.impl());
        blas::axpy(trans, SMD(1), BhA, C.impl());
    }
}

//
//  Case 3: trans==Trans, alpha
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
     const MatrixClosureR2KC3_<SV, MA1, MB1, MB2, MA2> &aAhB_aBhA,
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    const auto &aAhB = aAhB_aBhA.left();
    const auto &aBhA = aAhB_aBhA.right();

    const auto &a1 = aAhB.left().left().value();
    const auto &A1 = aAhB.left().right().right().right();
    const auto &B1 = aAhB.right();

    const auto &a2 = aBhA.left().left().value();
    const auto &B2 = aBhA.left().right().right().right();
    const auto &A2 = aBhA.right();

    if (identical(A1, A2)
     && identical(B1, B2)
     && (a1==cxxblas::conjugate(a2)))
    {
        trans = Transpose(trans^ConjTrans);
        blas::r2k(trans, alpha*a1, A1, B1, ALPHA(1), C.impl());
        return;
    } else {
        std::cerr << "warning" << std::endl;
        blas::axpy(trans, alpha, aAhB, C.impl());
        blas::axpy(trans, alpha, aBhA, C.impl());
    }
}

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
     const MatrixClosureR2KC3_<SV, MA1, MB1, MB2, MA2> &aAhB_aBhA,
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    const auto &aAhB = aAhB_aBhA.left();
    const auto &aBhA = aAhB_aBhA.right();

    const auto &a1 = aAhB.left().left().value();
    const auto &A1 = aAhB.left().right().right().right();
    const auto &B1 = aAhB.right();

    const auto &a2 = aBhA.left().left().value();
    const auto &B2 = aBhA.left().right().right().right();
    const auto &A2 = aBhA.right();

    typedef typename MC::Impl::ElementType  ALPHA;

    if (identical(A1, A2)
     && identical(B1, B2)
     && (a1==cxxblas::conjugate(a2)))
    {
        trans = Transpose(trans^ConjTrans);
        blas::r2k(trans, a1, A1, B1, ALPHA(0), C.impl());
        return;
    } else {
        std::cerr << "warning" << std::endl;
        blas::copy(trans, aAhB, C.impl());
        blas::axpy(trans, ALPHA(1), aBhA, C.impl());
    }
}

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
     const MatrixClosureR2KC3<MD, SV, MA1, MB1, MB2, MA2> &bD_aAhB_aBhA,
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MD>::Remainder     RMD;
    typedef typename PruneScaling<MD>::ScalingType   SMD;

    const auto &bD   = bD_aAhB_aBhA.left().left();
    const auto &aAhB = bD_aAhB_aBhA.left().right();
    const auto &aBhA = bD_aAhB_aBhA.right();

    const SMD &b    = PruneScaling<MD>::getFactor(bD);
    const RMD &D    = PruneScaling<MD>::getRemainder(bD);

    const auto &a1 = aAhB.left().left().value();
    const auto &A1 = aAhB.left().right().right().right();
    const auto &B1 = aAhB.right();

    const auto &a2 = aBhA.left().left().value();
    const auto &B2 = aBhA.left().right().right().right();
    const auto &A2 = aBhA.right();

    typedef typename MC::Impl::ElementType  ALPHA;

    if (identical(A1, A2)
     && identical(B1, B2)
     && (a1==cxxblas::conjugate(a2)))
    {
        trans = Transpose(trans^ConjTrans);
        blas::r2k(trans, a1, A1, B1, b, C.impl());
        return;
    } else {
        blas::copy(trans, D, C.impl());
        blas::scal(b, C.impl());
        blas::axpy(trans, ALPHA(1), aAhB, C.impl());
        blas::axpy(trans, ALPHA(1), aBhA, C.impl());
    }
}

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_TWEAKS_R2K_TCC
