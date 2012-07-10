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

#ifndef FLENS_BLAS_CLOSURES_MV_TCC
#define FLENS_BLAS_CLOSURES_MV_TCC 1

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//== TriangularMatrix - Vector products ========================================
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
void
mv(Transpose trans,
   const ALPHA &alpha, const TriangularMatrix<MA> &_A, const Vector<VX> &_x,
   const BETA &beta, Vector<VY> &_y)
{
    using namespace DEBUGCLOSURE;

    const typename TriangularMatrix<MA>::Impl  &A = _A.impl();
    const typename Vector<VX>::Impl            &x = _x.impl();
    typename Vector<VY>::Impl                  &y = _y.impl();

//
//  If beta is not Zero we need a temporary
//
#   ifndef FLENS_DEBUG_CLOSURES
    ASSERT(beta==BETA(0));
#   else
    const bool tmpNeeded = (beta!=BETA(0));
    typename Vector<VY>::Impl  yTmp;
    if (tmpNeeded) {
        FLENS_BLASLOG_TMP_ADD(yTmp);
        yTmp = y;
    }
#   endif

    if (!identical(x, y)) {
        y = x;
    }
    mv(trans, A, y);

    if (alpha!=ALPHA(1)) {
        scal(alpha, y.impl());
    }

#   ifdef FLENS_DEBUG_CLOSURES
    if (tmpNeeded) {
        y += beta*yTmp;
        FLENS_BLASLOG_TMP_REMOVE(yTmp, y);
    }
#   endif
}

//== SymmetricMatrix - Vector products =========================================
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
void
mv(Transpose trans,
   const ALPHA &alpha, const SymmetricMatrix<MA> &A, const Vector<VX> &x,
   const BETA &beta, Vector<VY> &y)
{
#   ifdef FLENS_DEBUG_CLOSURES

    if (trans==NoTrans || trans==Trans) {
        mv(alpha, A.impl(), x.impl(), beta, y.impl());
    } else {
        typedef typename MA::Impl::ElementType        TA;
        typedef GeMatrix<FullStorage<TA, ColMajor> >  RMA;

        RMA _A;
        FLENS_BLASLOG_TMP_ADD(_A);

        copy(trans, A, _A);

        FLENS_BLASLOG_TMP_REMOVE(_A, A);
        mv(NoTrans, alpha, _A, x.impl(), beta, y.impl());
    }

#   else

    ASSERT(trans==NoTrans || trans==Trans);
    mv(alpha, A.impl(), x.impl(), beta, y.impl());

#   endif
}

//== HermitianMatrix - Vector products =========================================
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
void
mv(Transpose trans,
   const ALPHA &alpha, const HermitianMatrix<MA> &A, const Vector<VX> &x,
   const BETA &beta, Vector<VY> &y)
{
#   ifdef FLENS_DEBUG_CLOSURES

    if (trans==NoTrans || trans==ConjTrans) {
        mv(alpha, A.impl(), x.impl(), beta, y.impl());
    } else {
        typedef typename MA::Impl::ElementType        TA;
        typedef GeMatrix<FullStorage<TA, ColMajor> >  RMA;

        RMA _A;
        FLENS_BLASLOG_TMP_ADD(_A);

        copy(trans, A, _A);

        FLENS_BLASLOG_TMP_REMOVE(_A, A);
        mv(NoTrans, alpha, _A, x.impl(), beta, y.impl());
    }

#   else

    ASSERT(trans==NoTrans || trans==ConjTrans);
    mv(alpha, A.impl(), x.impl(), beta, y.impl());

#   endif
}


//== Matrix - Vector products ==================================================
//
//  This gets called if everything else fails
//
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
void
mv(Transpose trans,
   const ALPHA &alpha, const Matrix<MA> &A, const Vector<VX> &x,
   const BETA &beta, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_GEMV(trans, alpha, A, x, beta, y);

    typedef typename MA::Impl::ElementType        TA;
    typedef typename VX::Impl::ElementType        TX;

    typedef GeMatrix<FullStorage<TA, ColMajor> >  RMA;
    typedef DenseVector<Array<TX> >               RVX;

    FLENS_BLASLOG_TMP_TRON;
    RMA _A = A.impl();
    RVX _x = x.impl();
    FLENS_BLASLOG_TMP_TROFF;

    mv(trans, alpha, _A, _x, beta, y.impl());

    FLENS_BLASLOG_TMP_REMOVE(_A, A);
    FLENS_BLASLOG_TMP_REMOVE(_x, x);

    FLENS_BLASLOG_END;
}

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_MV_TCC

