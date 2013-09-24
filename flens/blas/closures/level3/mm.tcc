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

#ifndef FLENS_BLAS_CLOSURES_LEVEL3_MM_TCC
#define FLENS_BLAS_CLOSURES_LEVEL3_MM_TCC 1

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

//== GeneralMatrix - GeneralMatrix products ====================================
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(Transpose transA, Transpose transB, const ALPHA &alpha,
   const GeneralMatrix<MA> &A, const GeneralMatrix<MB> &B,
   const BETA &beta, Matrix<MC> &C)
{
    // You get here if you want to call a matrix-matrix product that was not
    // defined.  Or its not correctly included.
    // TODO: print the signature of the blas::mm() function that needs to be
    //       implemented.
    ASSERT(0);
}

//== TriangularMatrix - GeneralMatrix products =================================
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
typename RestrictTo<IsTriangularMatrix<MA>::value &&
                    IsGeneralMatrix<MB>::value &&
                    IsGeneralMatrix<MC>::value,
         void>::Type
trmm(Side side, Transpose transA, Transpose transB, const ALPHA &alpha,
     const MA &_A, const MB &_B, const BETA &beta, MC &C)
{
    using namespace DEBUGCLOSURE;

    typedef typename Result<MA>::Type  RMA;
    typedef typename Result<MB>::Type  RMB;

//
//  In non-closure debug mode we do not allow temporaries for A or B.
//
#   ifndef FLENS_DEBUG_CLOSURES
    static_assert(IsSame<RMA, typename Result<RMA>::Type>::value,
                  "temporary required");
    static_assert(IsSame<RMB, typename Result<RMB>::Type>::value,
                  "temporary required");
#   else
    typedef typename Result<MC>::Type  RMC;
#   endif

//
//  If _A or _B is a closure temporaries get created
//
    FLENS_BLASLOG_TMP_TRON;
    const RMA &A = _A;
    const RMB &B = _B;
    FLENS_BLASLOG_TMP_TROFF;

//
//  If beta!=0 or transB!=NoTrans or A and C share the same memopry we need
//  temporaries
//
#   ifndef FLENS_DEBUG_CLOSURES
    ASSERT(beta==BETA(0));
    ASSERT(transB==NoTrans);
    ASSERT(!identical(A, C));
#   else
    if (transB!=NoTrans) {
//
//      apply op(B) and recall trmm
//
        typename RMB::NoView _B;
        FLENS_BLASLOG_TMP_ADD(_B);
        copy(transB, B, _B);
        trmm(side, transA, NoTrans, alpha, A, _B, beta, C);
        FLENS_BLASLOG_TMP_REMOVE(_B, B);
        if (!IsSame<RMA, typename Result<RMA>::Type>::value) {
            FLENS_BLASLOG_TMP_REMOVE(A, _A);
        }
        if (!IsSame<RMB, typename Result<RMB>::Type>::value) {
            FLENS_BLASLOG_TMP_REMOVE(B, _B);
        }
        return;
    }
    if (identical(A,C)) {
        FLENS_BLASLOG_IDENTICAL(A, C);
        typename RMC::NoView _C;
        FLENS_BLASLOG_TMP_ADD(_C);

        trmm(side, transA, NoTrans, alpha, A, B, beta, _C);
        C = _C;

        FLENS_BLASLOG_TMP_REMOVE(_C, C);
        return;
    }
    typename RMC::NoView tmpC;
    if (beta!=BETA(0)) {
        FLENS_BLASLOG_TMP_ADD(tmpC);
        tmpC = C;
    }
#   endif

//
//  trmm can only compute C = A*C or C = C*A.  So if B and C are not identical
//  we need to copy
//
    if (!identical(C, B)) {
        C = B;
    }
    mm(side, transA, alpha, A, C);

#   ifdef FLENS_DEBUG_CLOSURES
    if (beta!=BETA(0)) {
        C += beta*tmpC;
        FLENS_BLASLOG_TMP_REMOVE(tmpC, C);
    }
    if (!IsSame<RMA, typename Result<RMA>::Type>::value) {
        FLENS_BLASLOG_TMP_REMOVE(A, _A);
    }
    if (!IsSame<RMB, typename Result<RMB>::Type>::value) {
        FLENS_BLASLOG_TMP_REMOVE(B, _B);
    }
#   endif
}

template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(Transpose transA, Transpose transB, const ALPHA &alpha,
   const TriangularMatrix<MA> &A, const GeneralMatrix<MB> &B,
   const BETA &beta, Matrix<MC> &C)
{
    trmm(Left, transA, transB, alpha, A.impl(), B.impl(), beta, C.impl());
}

template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(Transpose transA, Transpose transB, const ALPHA &alpha,
   const GeneralMatrix<MA> &A, const TriangularMatrix<MB> &B,
   const BETA &beta, Matrix<MC> &C)
{
    trmm(Right, transB, transA, alpha, B.impl(), A.impl(), beta, C.impl());
}

//== SymmetricMatrix - GeneralMatrix products ==================================
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
symm(Side side, Transpose transB, const ALPHA &alpha,
     const SymmetricMatrix<MA> &_A, const GeneralMatrix<MB> &_B,
     const BETA &beta, Matrix<MC> &C)
{
    using namespace DEBUGCLOSURE;

    typedef typename Result<typename MA::Impl>::Type  RMA;
    typedef typename Result<typename MB::Impl>::Type  RMB;

//
//  In non-closure debug mode we do not allow temporaries for A or B.
//
#   ifndef FLENS_DEBUG_CLOSURES
    static_assert(IsSame<RMA, typename Result<RMA>::Type>::value,
                  "temporary required");
    static_assert(IsSame<RMB, typename Result<RMB>::Type>::value,
                  "temporary required");
    ASSERT(transB==NoTrans);
#   else
    typedef typename Result<typename MC::Impl>::Type  RMC;
#   endif

//
//  If _A or _B is a closure temporaries get created
//
    FLENS_BLASLOG_TMP_TRON;
    const RMA &A = _A.impl();
    const RMB &B = _B.impl();
    FLENS_BLASLOG_TMP_TROFF;

//
//  call (sy)mm
//
#   ifndef FLENS_DEBUG_CLOSURES
    mm(side, alpha, A.impl(), B.impl(), beta, C.impl());
#   else
//
//  if transB is not NoTrans we need another temporary
//
    if (transB==NoTrans) {
        mm(side, alpha, A, B, beta, C.impl());
    } else {
        typename RMB::NoView _B;
        FLENS_BLASLOG_TMP_ADD(_B);
        copy(transB, B, _B);
        mm(side, alpha, A, _B, beta, C.impl());
        FLENS_BLASLOG_TMP_REMOVE(_B, B);
    }
#   endif

#   ifdef FLENS_DEBUG_CLOSURES
    if (!IsSame<RMA, typename Result<RMA>::Type>::value) {
        FLENS_BLASLOG_TMP_REMOVE(A, _A);
    }
    if (!IsSame<RMB, typename Result<RMB>::Type>::value) {
        FLENS_BLASLOG_TMP_REMOVE(B, _B);
    }
#   endif
}

template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(Transpose transA, Transpose transB, const ALPHA &alpha,
   const SymmetricMatrix<MA> &A, const GeneralMatrix<MB> &B,
   const BETA &beta, Matrix<MC> &C)
{
    ASSERT(transA==NoTrans || transA==Trans);
    symm(Left, transB, alpha, A.impl(), B.impl(), beta, C.impl());
}

template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(Transpose transA, Transpose transB, const ALPHA &alpha,
   const GeneralMatrix<MA> &A, const SymmetricMatrix<MB> &B,
   const BETA &beta, Matrix<MC> &C)
{
    ASSERT(transA==NoTrans || transA==Trans);
    symm(Right, transB, alpha, B.impl(), A.impl(), beta, C.impl());
}

//== HermitianMatrix - GeneralMatrix products ==================================


template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
hemm(Side side, Transpose transB, const ALPHA &alpha,
     const HermitianMatrix<MA> &_A, const GeneralMatrix<MB> &_B,
     const BETA &beta, Matrix<MC> &C)
{
    using namespace DEBUGCLOSURE;

    typedef typename Result<typename MA::Impl>::Type  RMA;
    typedef typename Result<typename MB::Impl>::Type  RMB;

//
//  In non-closure debug mode we do not allow temporaries for A or B.
//
#   ifndef FLENS_DEBUG_CLOSURES
    static_assert(IsSame<RMA, typename Result<RMA>::Type>::value,
                  "temporary required");
    static_assert(IsSame<RMB, typename Result<RMB>::Type>::value,
                  "temporary required");
    ASSERT(transB==NoTrans);
#   else
    typedef typename Result<typename MC::Impl>::Type  RMC;
#   endif

//
//  If _A or _B is a closure temporaries get created
//
    FLENS_BLASLOG_TMP_TRON;
    const RMA &A = _A.impl();
    const RMB &B = _B.impl();
    FLENS_BLASLOG_TMP_TROFF;

//
//  call (he)mm
//
#   ifndef FLENS_DEBUG_CLOSURES
    mm(side, alpha, A.impl(), B.impl(), beta, C.impl());
#   else
//
//  if transB is not NoTrans we need another temporary
//
    if (transB==NoTrans) {
        mm(side, alpha, A, B, beta, C.impl());
    } else {
        typename RMB::NoView _B;
        FLENS_BLASLOG_TMP_ADD(_B);
        copy(transB, B, _B);
        mm(side, alpha, A, _B, beta, C.impl());
        FLENS_BLASLOG_TMP_REMOVE(_B, B);
    }
#   endif

#   ifdef FLENS_DEBUG_CLOSURES
    if (!IsSame<RMA, typename Result<RMA>::Type>::value) {
        FLENS_BLASLOG_TMP_REMOVE(A, _A);
    }
    if (!IsSame<RMB, typename Result<RMB>::Type>::value) {
        FLENS_BLASLOG_TMP_REMOVE(B, _B);
    }
#   endif
}


template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(Transpose transA, Transpose transB, const ALPHA &alpha,
   const HermitianMatrix<MA> &A, const GeneralMatrix<MB> &B,
   const BETA &beta, Matrix<MC> &C)
{
    ASSERT(transA==NoTrans || transA==Trans);
    hemm(Left, transB, alpha, A.impl(), B.impl(), beta, C.impl());
}

template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(Transpose transA, Transpose transB, const ALPHA &alpha,
   const GeneralMatrix<MA> &A, const HermitianMatrix<MB> &B,
   const BETA &beta, Matrix<MC> &C)
{
    ASSERT(transA==NoTrans || transA==Trans);
    hemm(Right, transB, alpha, B.impl(), A.impl(), beta, C.impl());
}






//== Matrix - Matrix products ==================================================
//
//  This gets called if everything else fails
//
#ifdef FLENS_DEBUG_CLOSURES

template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
void
mm(Transpose transA, Transpose transB, const ALPHA &alpha,
   const Matrix<MA> &A, const Matrix<MB> &B,
   const BETA &beta, Matrix<MC> &C)
{
    FLENS_BLASLOG_BEGIN_GEMM(transA, transB, alpha, A.impl(),
                             B.impl(), beta, C.impl());
//
//  We create temporaries of type GeMatrix for all matrices on the right
//  hand side.  If A and B can be converted to GeMatrix types this at
//  least does the desired compuation.
//
    typedef typename MA::Impl::ElementType        TA;
    typedef typename MB::Impl::ElementType        TB;

    typedef GeMatrix<FullStorage<TA, ColMajor> >  RMA;
    typedef GeMatrix<FullStorage<TB, ColMajor> >  RMB;

    FLENS_BLASLOG_TMP_TRON;
    const RMA &_A = A.impl();
    const RMB &_B = B.impl();
    FLENS_BLASLOG_TMP_TROFF;

    mm(transA, transB, alpha, _A, _B, beta, C.impl());

    FLENS_BLASLOG_TMP_REMOVE(_A, A);
    FLENS_BLASLOG_TMP_REMOVE(_B, B);

    FLENS_BLASLOG_END;
}

#endif

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_LEVEL3_MM_TCC

