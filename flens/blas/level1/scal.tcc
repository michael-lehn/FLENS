/*
 *   Copyright (c) 2009, Michael Lehn
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

#ifndef FLENS_BLAS_LEVEL1_SCAL_TCC
#define FLENS_BLAS_LEVEL1_SCAL_TCC 1

#include <cxxblas/cxxblas.h>
#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/closures/closures.h>
#include <flens/blas/level1/level1.h>
#include <flens/typedefs.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//-- BLAS Level 1 extensions ---------------------------------------------------

//-- scal
template <typename ALPHA, typename VY>
typename RestrictTo<IsDenseVector<VY>::value,
         void>::Type
scal(const ALPHA &alpha, VY &&y)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, y);

#   ifdef HAVE_CXXBLAS_SCAL
    cxxblas::scal(y.length(), alpha, y.data(), y.stride());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- scal
template <typename ALPHA, typename VY>
typename RestrictTo<IsTinyVector<VY>::value,
         void>::Type
scal(const ALPHA &alpha, VY &&y)
{
    typedef typename RemoveRef<VY>::Type   VectorY;
    typedef typename VectorY::ElementType  TY;

    const int n    = VectorY::Engine::length;
    const int incY = VectorY::Engine::stride;

    cxxblas::scal<n, ALPHA, TY, incY>(alpha, y.data());
}


//-- BLAS Level 1 extensions ---------------------------------------------------

//== GeneralMatrix

//-- gbscal
template <typename ALPHA, typename MB>
typename RestrictTo<IsDiagMatrix<MB>::value,
         void>::Type
scal(const ALPHA &alpha, MB &&B)
{
    scal(alpha, B.diag());
}
//-- gbscal
template <typename ALPHA, typename MB>
typename RestrictTo<IsGbMatrix<MB>::value,
         void>::Type
scal(const ALPHA &alpha, MB &&B)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

#   ifdef HAVE_CXXBLAS_GBSCAL
    cxxblas::gbscal(B.order(), B.numRows(), B.numCols(),
                    B.numSubDiags(), B.numSuperDiags(),
                    alpha, B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- gescal
template <typename ALPHA, typename MB>
typename RestrictTo<IsGeMatrix<MB>::value,
         void>::Type
scal(const ALPHA &alpha, MB &&B)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

#   ifdef HAVE_CXXBLAS_GESCAL
    cxxblas::gescal(B.order(), B.numRows(), B.numCols(),
                    alpha, B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- gescal
template <typename ALPHA, typename MB>
typename RestrictTo<IsGeTinyMatrix<MB>::value,
         void>::Type
scal(const ALPHA &alpha, MB &&B)
{
    typedef typename RemoveRef<MB>::Type   MatrixB;
    typedef typename MatrixB::ElementType  TB;

    const int m   = MatrixB::Engine::numRows;
    const int n   = MatrixB::Engine::numCols;
    const int ldB = MatrixB::Engine::leadingDimension;

    cxxblas::gescal<m,n,ALPHA,TB,ldB>(alpha, B.data());
}

//== HermitianMatrix

//-- hbscal
template <typename ALPHA, typename MB>
typename RestrictTo<IsHbMatrix<MB>::value,
         void>::Type
scal(const ALPHA &alpha, MB &&B)
{
    ASSERT(cxxblas::imag(alpha)==0);

    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

    typedef typename HbMatrix<MB>::IndexType  IndexType;

    const IndexType numSubDiags   = (B.upLo()==Upper) ? 0 : B.numOffDiags();
    const IndexType numSuperDiags = (B.upLo()==Upper) ? B.numOffDiags() : 0;

#   ifdef HAVE_CXXBLAS_GBSCAL
    cxxblas::gbscal(B.order(), B.dim(), B.dim(),
                    numSubDiags, numSuperDiags,
                    alpha, B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- hescal
template <typename ALPHA, typename MB>
typename RestrictTo<IsHeMatrix<MB>::value,
         void>::Type
scal(const ALPHA &alpha, MB &&B)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

#   ifdef HAVE_CXXBLAS_GESCAL
    cxxblas::hescal(B.order(), B.upLo(), B.dim(),
                    alpha, B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- hpscal
template <typename ALPHA, typename MB>
typename RestrictTo<IsHpMatrix<MB>::value,
         void>::Type
scal(const ALPHA &alpha, MB &&B)
{
    ASSERT(cxxblas::imag(alpha)==0);
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

#   ifdef HAVE_CXXBLAS_TPSCAL
    cxxblas::tpscal(B.order(), B.upLo(), B.diag(),
                    B.dim(),
                    alpha, B.data());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//== SymmetricMatrix

//-- sbscal
template <typename ALPHA, typename MB>
typename RestrictTo<IsSbMatrix<MB>::value,
         void>::Type
scal(const ALPHA &alpha, MB &&B)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

    typedef typename SbMatrix<MB>::IndexType  IndexType;

    const IndexType numSubDiags   = (B.upLo()==Upper) ? 0 : B.numOffDiags();
    const IndexType numSuperDiags = (B.upLo()==Upper) ? B.numOffDiags() : 0;

#   ifdef HAVE_CXXBLAS_GBSCAL
    cxxblas::gbscal(B.order(), B.dim(), B.dim(),
                    numSubDiags, numSuperDiags,
                    alpha, B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- spscal
template <typename ALPHA, typename MB>
typename RestrictTo<IsSpMatrix<MB>::value,
         void>::Type
scal(const ALPHA &alpha, MB &&B)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

#   ifdef HAVE_CXXBLAS_TPSCAL
    cxxblas::tpscal(B.order(), B.upLo(), B.diag(),
                    B.dim(),
                    alpha, B.data());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- syscal
template <typename ALPHA, typename MB>
typename RestrictTo<IsSyMatrix<MB>::value,
         void>::Type
scal(const ALPHA &alpha, MB &&B)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

#   ifdef HAVE_CXXBLAS_GESCAL
    cxxblas::syscal(B.order(), B.upLo(), B.dim(),
                    alpha, B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//== TriangularMatrix

//-- tbscal
template <typename ALPHA, typename MB>
typename RestrictTo<IsTbMatrix<MB>::value,
         void>::Type
scal(const ALPHA &alpha, MB &&B)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

#   ifdef HAVE_CXXBLAS_GBSCAL
    cxxblas::gbscal(B.order(), B.numRows(), B.numCols(),
                    B.numSubDiags(), B.numSuperDiags(),
                    alpha, B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- tpscal
template <typename ALPHA, typename MB>
typename RestrictTo<IsTpMatrix<MB>::value,
         void>::Type
scal(const ALPHA &alpha, MB &&B)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

#   ifdef HAVE_CXXBLAS_TPSCAL
    cxxblas::tpscal(B.order(), B.upLo(), B.diag(),
                    B.dim(),
                    alpha, B.data());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- trscal
template <typename ALPHA, typename MB>
typename RestrictTo<IsTrMatrix<MB>::value,
         void>::Type
scal(const ALPHA &alpha, MB &&B)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

    if (B.diag()==Unit) {
        if (B.upLo()==Upper) {
            scal(alpha, B.general().strictUpper());
        } else {
            scal(alpha, B.general().strictLower());
        }
        return;
    }

#   ifdef HAVE_CXXBLAS_TRSCAL
    typedef typename RemoveRef<MB>::Type   MatrixB;
    cxxblas::trscal(B.order(), B.upLo(), B.numRows(), B.numCols(),
                    alpha, B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_SCAL_TCC
