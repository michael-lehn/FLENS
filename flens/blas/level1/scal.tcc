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
#include <flens/auxiliary/macros.h>
#include <flens/storage/storage.h>
#include <flens/typedefs.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//-- scal (forwarding)
template <typename ALPHA, typename VY>
typename RestrictTo<IsSame<VY, typename VY::Impl>::value,
                    void>::Type
scal(const ALPHA &alpha, VY &&y)
{
    CHECKPOINT_ENTER;
    scal(alpha, y);
    CHECKPOINT_LEAVE;
}

//-- scal
template <typename ALPHA, typename VY>
void
scal(const ALPHA &alpha, DenseVector<VY> &y)
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

//-- gescal
template <typename ALPHA, typename MB>
void
scal(const ALPHA &alpha, GeMatrix<MB> &B)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

#   ifdef HAVE_CXXBLAS_GESCAL
    cxxblas::gescal(MB::order, B.numRows(), B.numCols(),
                    alpha, B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- gbscal
template <typename ALPHA, typename MB>
void
scal(const ALPHA &alpha, GbMatrix<MB> &B)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

#   ifdef HAVE_CXXBLAS_GBSCAL
    cxxblas::gbscal(MB::order, B.numRows(), B.numCols(),
                    B.numSubDiags(), B.numSuperDiags(),
                    alpha, B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- hbscal
template <typename ALPHA, typename MB>
void
scal(const ALPHA &alpha, HbMatrix<MB> &B)
{
    ASSERT(cxxblas::imag(alpha)==0);
    
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

    typedef typename HbMatrix<MB>::IndexType  IndexType;
    
    const IndexType numSubDiags   = (B.upLo()==Upper) ? 0 : B.numOffDiags();
    const IndexType numSuperDiags = (B.upLo()==Upper) ? B.numOffDiags() : 0;
    
#   ifdef HAVE_CXXBLAS_GBSCAL
    cxxblas::gbscal(MB::order, B.dim(), B.dim(),
                    numSubDiags, numSuperDiags,
                    alpha, B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- hpscal
template <typename ALPHA, typename MB>
void
scal(const ALPHA &alpha, HpMatrix<MB> &B)
{
    ASSERT(cxxblas::imag(alpha)==0);
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

#   ifdef HAVE_CXXBLAS_TPSCAL
    cxxblas::tpscal(MB::order, B.upLo(), B.diag(),
                    B.dim(), 
                    alpha, B.data());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- sbscal
template <typename ALPHA, typename MB>
void
scal(const ALPHA &alpha, SbMatrix<MB> &B)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

    typedef typename SbMatrix<MB>::IndexType  IndexType;
    
    const IndexType numSubDiags   = (B.upLo()==Upper) ? 0 : B.numOffDiags();
    const IndexType numSuperDiags = (B.upLo()==Upper) ? B.numOffDiags() : 0;
    
#   ifdef HAVE_CXXBLAS_GBSCAL
    cxxblas::gbscal(MB::order, B.dim(), B.dim(),
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
void
scal(const ALPHA &alpha, SpMatrix<MB> &B)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

#   ifdef HAVE_CXXBLAS_TPSCAL
    cxxblas::tpscal(MB::order, B.upLo(), B.diag(),
                    B.dim(), 
                    alpha, B.data());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- tbscal
template <typename ALPHA, typename MB>
void
scal(const ALPHA &alpha, TbMatrix<MB> &B)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

#   ifdef HAVE_CXXBLAS_GBSCAL
    cxxblas::gbscal(MB::order, B.numRows(), B.numCols(),
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
void
scal(const ALPHA &alpha, TpMatrix<MB> &B)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_SCAL(alpha, B);

#   ifdef HAVE_CXXBLAS_TPSCAL
    cxxblas::tpscal(MB::order, B.upLo(), B.diag(),
                    B.dim(), 
                    alpha, B.data());
#   else
    ASSERT(0);
#   endif

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_SCAL_TCC
