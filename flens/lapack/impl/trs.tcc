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
       SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --

       SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
       SUBROUTINE ZTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 *
 */

#ifndef FLENS_LAPACK_IMPL_TRS_TCC
#define FLENS_LAPACK_IMPL_TRS_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- (ge)trs [real and complex variant] ----------------------------------------

template <typename MA, typename VP, typename MB>
void
trs_impl(Transpose trans, const GeMatrix<MA> &A, const DenseVector<VP> &piv,
         GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType    IndexType;
    typedef typename GeMatrix<MA>::ElementType  T;

    const IndexType n       = A.numCols();
    const IndexType nRhs    = B.numCols();

    const T  One(1);
//
//  Quick return if possible
//
    if ((n==0) || (nRhs==0)) {
        return;
    }

    if ((trans==NoTrans) || (trans==Conj)) {
//
//      Solve A * X = B.
//
//      Apply row interchanges to the right hand sides.
//
        laswp(B, piv);
//
//      Solve L*X = B, overwriting B with X.
//
        blas::sm(Left, trans, One, A.lowerUnit(), B);
//
//      Solve U*X = B, overwriting B with X.
//
        blas::sm(Left, trans, One, A.upper(), B);
    } else {
//
//      Solve A' * X = B.
//
//      Solve U'*X = B, overwriting B with X.
//
        blas::sm(Left, trans, One, A.upper(), B);
//
//      Solve L'*X = B, overwriting B with X.
//
        blas::sm(Left, trans, One, A.lowerUnit(), B);
//
//      Apply row interchanges to the solution vectors.
//
        laswp(B, piv.reverse());
    }
}

//-- (tr)trs [real and complex variant] ----------------------------------------

template <typename MA, typename MB>
typename TrMatrix<MA>::IndexType
trs_impl(Transpose trans, const TrMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename TrMatrix<MA>::IndexType    IndexType;
    typedef typename TrMatrix<MA>::ElementType  T;

    const IndexType n       = A.dim();

    const T  Zero(0), One(1);

    IndexType info = 0;
//
//  Quick return if possible
//
    if (n==0) {
        return info;
    }
//
//  Check for singularity.
//
    if (A.diag()!=Unit) {
        for (info=1; info<=n; ++info) {
            if (A(info,info)==Zero) {
                return info;
            }
        }
    }
    info = 0;
//
//  Solve A * x = b  or  A**T * x = b.
//
    blas::sm(Left, trans, One, A, B);

    return info;
}

} // namespace generic


//== interface for external lapack =============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (ge)trs [real and complex variant] ----------------------------------------

template <typename MA, typename VP, typename MB>
void
trs_impl(Transpose trans, const GeMatrix<MA> &A, const DenseVector<VP> &piv,
         GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info;
    info = cxxlapack::getrs<IndexType>(getF77Char(trans),
                                       A.numRows(),
                                       B.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       piv.data(),
                                       B.data(),
                                       B.leadingDimension());
    ASSERT(info==0);
}

//-- (tr)trs [real and complex variant] ----------------------------------------

template <typename MA, typename MB>
typename TrMatrix<MA>::IndexType
trs_impl(Transpose trans, const TrMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename TrMatrix<MA>::IndexType  IndexType;

    IndexType info;
    info = cxxlapack::trtrs<IndexType>(getF77Char(A.upLo()),
                                       getF77Char(trans),
                                       getF77Char(A.diag()),
                                       A.dim(),
                                       B.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       B.data(),
                                       B.leadingDimension());
    ASSERT(info>=0);
    return info;
}

//-- (he)trs [complex variant] ----------------------------------------

template <typename MA, typename VP, typename MB>
void
trs_impl(const HeMatrix<MA> &A, const DenseVector<VP> &piv,
         GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info;
    info = cxxlapack::hetrs<IndexType>(getF77Char(A.upLo()),
                                       A.dim(),
                                       B.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       piv.data(),
                                       B.data(),
                                       B.leadingDimension());
    ASSERT(info==0);
}

//-- (sy)trs [real variant] ----------------------------------------

template <typename MA, typename VP, typename MB>
void
trs_impl(const SyMatrix<MA> &A, const DenseVector<VP> &piv,
         GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info;
    info = cxxlapack::sytrs<IndexType>(getF77Char(A.upLo()),
                                       A.dim(),
                                       B.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       piv.data(),
                                       B.data(),
                                       B.leadingDimension());
    ASSERT(info==0);
}

//-- (gb)trs [real and complex variant] ----------------------------------------

template <typename MA, typename VP, typename MB>
void
trs_impl(Transpose trans, const GbMatrix<MA> &A, const DenseVector<VP> &piv,
         GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info;
    info = cxxlapack::gbtrs<IndexType>(getF77Char(trans),
                                       A.numRows(),
                                       A.numSubDiags(),
                                       A.numSuperDiags()-A.numSubDiags(),
                                       B.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       piv.data(),
                                       B.data(),
                                       B.leadingDimension());
    ASSERT(info==0);
}

//-- (tb)trs [real and complex variant] ----------------------------------------

template <typename MA, typename VP, typename MB>
void
trs_impl(Transpose trans, const TbMatrix<MA> &A, const DenseVector<VP> &piv,
         GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info;
    info = cxxlapack::tbtrs<IndexType>(getF77Char(A.upLo()),
                                       getF77Char(trans),
                                       A.dim(),
                                       A.numOffDiags(),
                                       B.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       B.data(),
                                       B.leadingDimension());
    ASSERT(info==0);
}

//-- (hp)trs [complex variant] ----------------------------------------

template <typename MA, typename VP, typename MB>
void
trs_impl(const HpMatrix<MA> &A, const DenseVector<VP> &piv,
         GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info;
    info = cxxlapack::hptrs<IndexType>(getF77Char(A.upLo()),
                                       A.dim(),
                                       B.numCols(),
                                       A.data(),
                                       piv.data(),
                                       B.data(),
                                       B.leadingDimension());
    ASSERT(info==0);
}

//-- (sp)trs [real variant] ----------------------------------------

template <typename MA, typename VP, typename MB>
void
trs_impl(const SpMatrix<MA> &A, const DenseVector<VP> &piv,
         GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info;
    info = cxxlapack::sptrs<IndexType>(getF77Char(A.upLo()),
                                       A.dim(),
                                       B.numCols(),
                                       A.data(),
                                       piv.data(),
                                       B.data(),
                                       B.leadingDimension());
    ASSERT(info==0);
}

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

//-- (tp)trs [real and complex variant] ----------------------------------------

template <typename MA, typename MB>
typename TpMatrix<MA>::IndexType
trs_impl(Transpose trans, const TpMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename TpMatrix<MA>::IndexType  IndexType;

    IndexType info;
    info = cxxlapack::tbtrs<IndexType>(getF77Char(A.upLo()),
                                       getF77Char(trans),
                                       getF77Char(A.diag()),
                                       A.dim(),
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

//-- (ge)trs [real and complex variant] ----------------------------------------

template <typename MA, typename VPIV, typename MB>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
trs(Transpose trans, const MA &A, const VPIV &piv, MB &&B)
{
    LAPACK_DEBUG_OUT("(ge)trs [real/complex]");
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type     MatrixA;
    typedef typename MatrixA::IndexType      IndexType;
    typedef typename RemoveRef<MB>::Type     MatrixB;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());

    const IndexType n = A.numRows();

    ASSERT(piv.firstIndex()==1);
    ASSERT(piv.length()==n);

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==n);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename MatrixB::NoView  B_org   = B;
#   endif    
    
//
//  Call implementation
//
    LAPACK_SELECT::trs_impl(trans, A, piv, B);
//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixB::NoView  B_generic   = B;

    B   = B_org;

    external::trs_impl(trans, A, piv, B);

    bool failed = false;
    if (! isIdentical(B_generic, B, "B_generic", "B")) {
        std::cerr << "CXXLAPACK: B_generic = " << B_generic << std::endl;
        std::cerr << "F77LAPACK: B = " << B << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    } else {
        // std::cerr << "passed: (ge)trs.tcc" << std::endl;
    }
#   endif
}

//-- (ge/gb)trs [variant if rhs is vector] ----------------------------------------

template <typename MA, typename VPIV, typename VB>
typename RestrictTo< (IsGeMatrix<MA>::value ||
                      IsGbMatrix<MA>::value )
                 && IsIntegerDenseVector<VPIV>::value
                 && IsDenseVector<VB>::value,
         void>::Type
trs(Transpose trans, const MA &A, const VPIV &piv, VB &&b)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type     MatrixA;
    typedef typename RemoveRef<VB>::Type     VectorB;

//
//  Create matrix view from vector b and call above variant
//
    typedef typename VectorB::ElementType  ElementType;
    typedef typename VectorB::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = MatrixA::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b, n);

    trs(trans, A, piv, B);
}

//-- (tr)trs [real and complex variant] ----------------------------------------

template <typename MA, typename MB>
typename RestrictTo<IsTrMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
trs(Transpose trans, const MA &A, MB &&B)
{
    LAPACK_DEBUG_OUT("(tr)trs [real/complex]");
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename RemoveRef<MB>::Type    MatrixB;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    const IndexType n = A.dim();

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==n);
#   endif

#   ifdef CHECK_CXXLAPACK
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

    IndexType _info = external::trs_impl(trans, A, B);

    bool failed = false;
    if (! isIdentical(B_generic, B, "B_generic", "B")) {
        std::cerr << "CXXLAPACK: B_generic = " << B_generic << std::endl;
        std::cerr << "F77LAPACK: B = " << B << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, "info", "_info")) {
        std::cerr << "CXXLAPACK: info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
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

//-- (tr/tb/tp)trs [variant if rhs is vector] ----------------------------------------

template <typename MA, typename VB>
typename RestrictTo< (IsTrMatrix<MA>::value ||
                      IsTbMatrix<MA>::value ||
                      IsTpMatrix<MA>::value)
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


//-- (he/sy/hp/sp)trs [variant if rhs is vector] ----------------------------------------

template <typename MA, typename VPIV, typename VB>
typename RestrictTo< (IsHeMatrix<MA>::value ||
                      IsSyMatrix<MA>::value ||
                      IsHpMatrix<MA>::value ||
                      IsSpMatrix<MA>::value )
                  && IsIntegerDenseVector<VPIV>::value
                  && IsDenseVector<VB>::value,
          void>::Type
trs(const MA &A, const VPIV &piv, VB &&b)
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

    trs(A, piv, B);
}

#ifdef USE_CXXLAPACK
//-- (he)trs [complex variant]

template <typename MA, typename VPIV, typename MB>
typename RestrictTo<IsHeMatrix<MA>::value
                && IsIntegerDenseVector<VPIV>::value
                && IsComplexGeMatrix<MB>::value,
        void>::Type
trs(const MA &A, const VPIV &piv, MB &&B)
{
    LAPACK_DEBUG_OUT("(he)trs [complex]");
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type     MatrixA;
    typedef typename MatrixA::IndexType      IndexType;
    typedef typename RemoveRef<MB>::Type     MatrixB;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());

    const IndexType n = A.numRows();

    ASSERT(piv.firstIndex()==1);
    ASSERT(piv.length()==n);

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==n);
#   endif

//
//  Call implementation
//
    external::trs_impl(A, piv, B);
}

//-- (sy)trs [real and complex variant]
template <typename MA, typename VPIV, typename MB>
typename RestrictTo<IsSyMatrix<MA>::value
                && IsIntegerDenseVector<VPIV>::value
                && IsGeMatrix<MB>::value,
        void>::Type
trs(const MA &A, const VPIV &piv, MB &&B)
{
    LAPACK_DEBUG_OUT("(sy)trs [real/complex]");
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type     MatrixA;
    typedef typename MatrixA::IndexType      IndexType;
    typedef typename RemoveRef<MB>::Type     MatrixB;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());

    const IndexType n = A.numRows();

    ASSERT(piv.firstIndex()==1);
    ASSERT(piv.length()==n);

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==n);
#   endif

//
//  Call implementation
//
    external::trs_impl(A, piv, B);
}

//-- (gb)trs [real and complex variant] ----------------------------------------

template <typename MA, typename VPIV, typename MB>
typename RestrictTo<IsGbMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
trs(Transpose trans, const MA &A, const VPIV &piv, MB &&B)
{
    LAPACK_DEBUG_OUT("(gb)trs [real/complex]");
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type     MatrixA;
    typedef typename MatrixA::IndexType      IndexType;
    typedef typename RemoveRef<MB>::Type     MatrixB;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());

    const IndexType n = A.numRows();

    ASSERT(piv.firstIndex()==1);
    ASSERT(piv.length()==n);

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==n);
#   endif

//
//  Call implementation
//
    external::trs_impl(trans, A, piv, B);

}

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
    typedef typename RemoveRef<MB>::Type    MatrixB;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    const IndexType n = A.dim();

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==n);
#   endif

//
//  Call implementation
//
    IndexType info = external::trs_impl(trans, A, B);

    return info;
}

//-- (hp)trs [complex variant]
template <typename MA, typename VPIV, typename MB>
typename RestrictTo<IsHpMatrix<MA>::value
		  && IsIntegerDenseVector<VPIV>::value
		  && IsComplexGeMatrix<MB>::value,
	  void>::Type
trs(const MA &A, const VPIV &piv, MB &&B)
{
    LAPACK_DEBUG_OUT("(hp)trs [complex]");
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type     MatrixA;
    typedef typename MatrixA::IndexType      IndexType;
    typedef typename RemoveRef<MB>::Type     MatrixB;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstIndex()==1);

    const IndexType n = A.dim();

    ASSERT(piv.firstIndex()==1);
    ASSERT(piv.length()==n);

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==n);
#   endif

//
//  Call implementation
//
    external::trs_impl(A, piv, B);
    
}

//-- (sp)trs [real and complex variant]
template <typename MA, typename VPIV, typename MB>
typename RestrictTo<IsSpMatrix<MA>::value
		  && IsIntegerDenseVector<VPIV>::value
		  && IsGeMatrix<MB>::value,
	  void>::Type
trs(const MA &A, const VPIV &piv, MB &&B)
{
    LAPACK_DEBUG_OUT("(sp)trs [real/complex]");
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type     MatrixA;
    typedef typename MatrixA::IndexType      IndexType;
    typedef typename RemoveRef<MB>::Type     MatrixB;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstIndex()==1);

    const IndexType n = A.dim();

    ASSERT(piv.firstIndex()==1);
    ASSERT(piv.length()==n);

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==n);
#   endif

//
//  Call implementation
//
    external::trs_impl(A, piv, B);

}

//-- (tp)trs [real and complex variant] ----------------------------------------

template <typename MA, typename MB>
typename RestrictTo<IsTpMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
trs(Transpose trans, const MA &A, MB &&B)
{
    LAPACK_DEBUG_OUT("(tp)trs [real/complex]");
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename RemoveRef<MB>::Type    MatrixB;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    const IndexType n = A.dim();

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==n);
#   endif

//
//  Call implementation
//
    IndexType info = external::trs_impl(trans, A, B);

    return info;
}


#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_TRS_TCC
