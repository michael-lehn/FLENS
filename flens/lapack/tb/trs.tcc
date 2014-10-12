/*
 *   Copyright (c) 2012, Michael Lehn, Klaus Pototzky
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
       SUBROUTINE DTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,
      $                   LDB, INFO )
       SUBROUTINE ZTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,
      $                   LDB, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 *
 */

#ifndef FLENS_LAPACK_TB_TRS_TCC
#define FLENS_LAPACK_TB_TRS_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- (tb)trs [real and complex variant] ----------------------------------------

template <typename MA, typename MB>
typename TbMatrix<MA>::IndexType
trs_impl(Transpose trans, const TbMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename TbMatrix<MA>::IndexType       IndexType;
    typedef typename TbMatrix<MA>::ElementType     ElementType;

    const ElementType zero(0);

    const Underscore<IndexType> _;

//
//     Check for singularity.
//
    if( A.diag()==NonUnit ) {
        for (IndexType i = 1; i<=A.dim(); ++i) {
            if( A(i, i)==zero ) {
                return i;
            }
        }
    }
//
//     Solve A * X = B,  A**T * X = B,  or  A**H * X = B.
//
    for (IndexType j=1; j<=B.numCols(); ++j) {
        blas::sv(trans, A, B(_,j));
    }
    return 0;
}

} // namespace generic

//== interface for external lapack =============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (tb)trs [real and complex variant] ----------------------------------------

template <typename MA, typename MB>
typename TbMatrix<MA>::IndexType
trs_impl(Transpose trans, const TbMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename TbMatrix<MA>::IndexType  IndexType;

    IndexType info;
    info = cxxlapack::tbtrs<IndexType>(getF77Char(A.upLo()),
                                       getF77Char(trans),
                                       getF77Char(A.diag()),
                                       A.dim(),
                                       A.numOffDiags(),
                                       B.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       B.data(),
                                       B.leadingDimension());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- (tb)trs [real and complex variant] ----------------------------------------

template <typename MA, typename MB>
typename RestrictTo<IsTbMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
trs(Transpose trans, const MA &A, MB &&B)
{
    LAPACK_DEBUG_OUT("(tb)trs [real/complex]");
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstIndex()==1);

    const IndexType n = A.dim();

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==n);
#   endif

#   ifdef CHECK_CXXLAPACK

    typedef typename RemoveRef<MB>::Type    MatrixB;
//
//  Make copies of output arguments
//
    typename MatrixB::NoView  B_org   = B;
#   endif
//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::trs_impl(trans, A, B);
//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixB::NoView  B_generic   = B;

    B   = B_org;

    IndexType info_ = external::trs_impl(trans, A, B);

    bool failed = false;
    if (! isIdentical(B_generic, B, "B_generic", "B")) {
        std::cerr << "CXXLAPACK: B_generic = " << B_generic << std::endl;
        std::cerr << "F77LAPACK: B = " << B << std::endl;
        failed = true;
    }

    if (! isIdentical(info, info_, "info", "info_")) {
        std::cerr << "CXXLAPACK: info = " << info << std::endl;
        std::cerr << "F77LAPACK: info_ = " << info_ << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    } else {
        // std::cerr << "passed: (tr)trs.tcc" << std::endl;
    }
#   endif

    return info;
}

//-- (tb)trs [variant if rhs is vector] ----------------------------------------

template <typename MA, typename VB>
typename RestrictTo<IsTbMatrix<MA>::value
                 && IsDenseVector<VB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
trs(Transpose trans, const MA &A, VB &&b)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VB>::Type    VectorB;

//
//  Create matrix view from vector b and call above variant
//
    typedef typename VectorB::ElementType  ElementType;
    typedef typename VectorB::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = MatrixA::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b, n);

    return trs(trans, A, B);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_TB_TRS_TCC
