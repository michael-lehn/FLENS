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

#ifndef FLENS_BLAS_CLOSURES_TWEAKS_RK_TCC
#define FLENS_BLAS_CLOSURES_TWEAKS_RK_TCC 1

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MA1>::Remainder    RMA;
    typedef typename PruneScaling<MA1>::ScalingType  SMA;

    const RMA &A      = PruneScaling<MA1>::getRemainder(aAAt.left());
    const SMA &alpha_ = PruneScaling<MA1>::getFactor(aAAt.left());

#   ifndef NDEBUG
    const auto &A_     = aAAt.right().right();
    ASSERT(identical(A, A_));
#   endif

    blas::rk(trans, alpha*alpha_, A, ALPHA(1), C.impl());
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MA1>::Remainder    RMA;
    typedef typename PruneScaling<MA1>::ScalingType  SMA;

    const RMA &A     = PruneScaling<MA1>::getRemainder(aAAt.left());
    const SMA &alpha = PruneScaling<MA1>::getFactor(aAAt.left());

#   ifndef NDEBUG
    const auto &A_   = aAAt.right().right();
    ASSERT(identical(A, A_));
#   endif

    typedef typename MC::Impl::ElementType  ALPHA;

    blas::rk(trans, alpha, A, ALPHA(0), C.impl());
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MB>::Remainder    RMB;
    typedef typename PruneScaling<MB>::ScalingType  SMB;

    typedef typename PruneScaling<MA1>::Remainder    RMA;
    typedef typename PruneScaling<MA1>::ScalingType  SMA;

    const RMB &B    = PruneScaling<MB>::getRemainder(bB_aAAt.left());
    const SMB &beta = PruneScaling<MB>::getFactor(bB_aAAt.left());

    const RMA &A     = PruneScaling<MA1>::getRemainder(bB_aAAt.right().left());
    const SMA &alpha = PruneScaling<MA1>::getFactor(bB_aAAt.right().left());

#   ifndef NDEBUG
    const auto &A_     = bB_aAAt.right().right().right();
    ASSERT(identical(A, A_));
#   endif

    if (identical(B, C.impl())) {
        blas::rk(trans, alpha, A, beta, C.impl());
    } else {
        blas::copy(trans, B, C.impl());
        blas::scal(beta, C.impl());
        blas::axpy(trans, SMA(1), bB_aAAt.right(), C.impl());
    }
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    const auto &A      = AtA.left().left();

#   ifndef NDEBUG
    const auto &A_     = AtA.right();
    ASSERT(identical(A, A_));
#   endif

    trans = Transpose(trans^Trans);
    blas::rk(trans, alpha, A, ALPHA(1), C.impl());
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    const auto &A      = AtA.left().left();

#   ifndef NDEBUG
    const auto &A_     = AtA.right();
    ASSERT(identical(A, A_));
#   endif

    typedef typename MC::Impl::ElementType  ALPHA;

    trans = Transpose(trans^Trans);
    blas::rk(trans, ALPHA(1), A, ALPHA(0), C.impl());
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MB>::Remainder    RMB;
    typedef typename PruneScaling<MB>::ScalingType  SMB;

    const RMB &B    = PruneScaling<MB>::getRemainder(bB_AtA.left());
    const SMB &beta = PruneScaling<MB>::getFactor(bB_AtA.left());

    const auto &A   = bB_AtA.right().left().left();

#   ifndef NDEBUG
    const auto &A_  = bB_AtA.right().right();
    ASSERT(identical(A, A_));
#   endif

    trans = Transpose(trans^Trans);
    if (identical(B, C.impl())) {
        blas::rk(trans, SMB(1), A, beta, C.impl());
    } else {
        blas::copy(trans, B, C.impl());
        blas::scal(beta, C.impl());
        blas::axpy(trans, SMB(1), bB_AtA.right(), C.impl());
    }
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;


    const auto &alpha_ = aAtA.left().left().value();
    const auto &A      = aAtA.left().right().right();

#   ifndef NDEBUG
    const auto &A_    = aAtA.right();
    ASSERT(identical(A, A_));
#   endif

    trans = Transpose(trans^Trans);
    blas::rk(trans, alpha*alpha_, A, ALPHA(1), C.impl());
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;


    const auto &alpha = aAtA.left().left().value();
    const auto &A     = aAtA.left().right().right();

#   ifndef NDEBUG
    const auto &A_    = aAtA.right();
    ASSERT(identical(A, A_));
#   endif

    typedef typename MC::Impl::ElementType  ALPHA;

    trans = Transpose(trans^Trans);
    blas::rk(trans, alpha, A, ALPHA(0), C.impl());
}

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
     const MatrixClosureRKU3<MB, SV, MA1, MA2> &bB_aAtA,
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=ConjTrans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MB>::Remainder    RMB;
    typedef typename PruneScaling<MB>::ScalingType  SMB;

    const RMB &B      = PruneScaling<MB>::getRemainder(bB_aAtA.left());
    const SMB &beta   = PruneScaling<MB>::getFactor(bB_aAtA.left());

    const auto &alpha = bB_aAtA.right().left().left().value();
    const auto &A     = bB_aAtA.right().left().right().right();

#   ifndef NDEBUG
    const auto &A_     = bB_aAtA.right().right();
    ASSERT(identical(A, A_));
#   endif

    trans = Transpose(trans^Trans);
    if (identical(B, C.impl())) {
        blas::rk(trans, alpha, A, beta, C.impl());
    } else {
        blas::copy(trans, B, C.impl());
        blas::scal(beta, C.impl());
        blas::axpy(trans, alpha, bB_aAtA.right(), C.impl());
    }
}


//-- HermitianMatrix -----------------------------------------------------------
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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MA1>::Remainder    RMA;
    typedef typename PruneScaling<MA1>::ScalingType  SMA;

    const RMA &A      = PruneScaling<MA1>::getRemainder(aAAh.left());
    const SMA &alpha_ = PruneScaling<MA1>::getFactor(aAAh.left());

#   ifndef NDEBUG
    const auto &A_     = aAAh.right().right().right();
    ASSERT(identical(A, A_));
#   endif

    blas::rk(trans, alpha*alpha_, A, ALPHA(1), C.impl());
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MA1>::Remainder    RMA;
    typedef typename PruneScaling<MA1>::ScalingType  SMA;

    const RMA &A     = PruneScaling<MA1>::getRemainder(aAAh.left());
    const SMA &alpha = PruneScaling<MA1>::getFactor(aAAh.left());

#   ifndef NDEBUG
    const auto &A_     = aAAh.right().right().right();
    ASSERT(identical(A, A_));
#   endif

    typedef typename MC::Impl::ElementType  ALPHA;

    blas::rk(trans, alpha, A, ALPHA(0), C.impl());
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MB>::Remainder    RMB;
    typedef typename PruneScaling<MB>::ScalingType  SMB;

    typedef typename PruneScaling<MA1>::Remainder    RMA;
    typedef typename PruneScaling<MA1>::ScalingType  SMA;

    const auto &bB   = bB_aAAh.left();
    const auto &aAAh = bB_aAAh.right();

    const RMB &B     = PruneScaling<MB>::getRemainder(bB);
    const SMB &beta  = PruneScaling<MB>::getFactor(bB);

    const SMA &alpha = PruneScaling<MA1>::getFactor(aAAh.left());
    const RMA &A     = PruneScaling<MA1>::getRemainder(aAAh.left());

#   ifndef NDEBUG
    const auto &A_     = aAAh.right().right().right();
    ASSERT(identical(A, A_));
#   endif

    if (identical(B, C.impl())) {
        blas::rk(trans, alpha, A, beta, C.impl());
    } else {
        blas::copy(trans, B, C.impl());
        blas::scal(beta, C.impl());
        blas::axpy(trans, SMA(1), aAAh, C.impl());
    }
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    const auto &A      = AhA.left().left().left();

#   ifndef NDEBUG
    const auto &A_     = AhA.right();
    ASSERT(identical(A, A_));
#   endif

    trans = Transpose(trans^ConjTrans);
    blas::rk(trans, alpha, A, ALPHA(1), C.impl());
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    const auto &A      = AhA.left().left().left();

#   ifndef NDEBUG
    const auto &A_     = AhA.right();
    ASSERT(identical(A, A_));
#   endif

    typedef typename MC::Impl::ElementType  ALPHA;

    trans = Transpose(trans^ConjTrans);
    blas::rk(trans, ALPHA(1), A, ALPHA(0), C.impl());
}

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
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MB>::Remainder    RMB;
    typedef typename PruneScaling<MB>::ScalingType  SMB;

    const auto &bB  = bB_AhA.left();
    const auto &AhA = bB_AhA.right();

    const RMB &B    = PruneScaling<MB>::getRemainder(bB);
    const SMB &beta = PruneScaling<MB>::getFactor(bB);

    const auto &A   = AhA.left().left().left();

#   ifndef NDEBUG
    const auto &A_  = AhA.right();
    ASSERT(identical(A, A_));
#   endif

    trans = Transpose(trans^ConjTrans);
    if (identical(B, C.impl())) {
        blas::rk(trans, SMB(1), A, beta, C.impl());
    } else {
        blas::copy(trans, B, C.impl());
        blas::scal(beta, C.impl());
        blas::axpy(trans, SMB(1), AhA, C.impl());
    }
}

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
     const MatrixClosureRKC3_<SV, MA1, MA2> &aAhA,
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;


    const auto &alpha_ = aAhA.left().left().value();
    const auto &A      = aAhA.left().right().right().right();

#   ifndef NDEBUG
    const auto &A_    = aAhA.right();
    ASSERT(identical(A, A_));
#   endif

    trans = Transpose(trans^ConjTrans);
    blas::rk(trans, alpha*alpha_, A, ALPHA(1), C.impl());
}

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
     const MatrixClosureRKC3_<SV, MA1, MA2> &aAhA,
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;


    const auto &alpha = aAhA.left().left().value();
    const auto &A     = aAhA.left().right().right().right();

#   ifndef NDEBUG
    const auto &A_    = aAhA.right();
    ASSERT(identical(A, A_));
#   endif

    typedef typename MC::Impl::ElementType  ALPHA;

    trans = Transpose(trans^ConjTrans);
    blas::rk(trans, alpha, A, ALPHA(0), C.impl());
}

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
     const MatrixClosureRKC3<MB, SV, MA1, MA2> &bB_aAhA,
     MC &C)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans!=Conj);
    ASSERT(trans!=Trans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<MB>::Remainder    RMB;
    typedef typename PruneScaling<MB>::ScalingType  SMB;

    const auto &bB   = bB_aAhA.left();
    const auto &aAhA = bB_aAhA.right();

    const RMB &B      = PruneScaling<MB>::getRemainder(bB);
    const SMB &beta   = PruneScaling<MB>::getFactor(bB);

    const auto &alpha = aAhA.left().left().value();
    const auto &A     = aAhA.left().right().right().right();

#   ifndef NDEBUG
    const auto &A_     = aAhA.right();
    ASSERT(identical(A, A_));
#   endif

    trans = Transpose(trans^ConjTrans);
    if (identical(B, C.impl())) {
        blas::rk(trans, alpha, A, beta, C.impl());
    } else {
        blas::copy(trans, B, C.impl());
        blas::scal(beta, C.impl());
        blas::axpy(trans, alpha, aAhA, C.impl());
    }
}

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_TWEAKS_RK_TCC
