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
       SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
     $                  LWORK, INFO )
       SUBROUTINE ZHESV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
     $                  LWORK, INFO )
       SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
       SUBROUTINE ZGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
       SUBROUTINE DSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
       SUBROUTINE ZSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
 *
 *  -- LAPACK driver routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_SV_TCC
#define FLENS_LAPACK_IMPL_SV_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- (ge)sv [real and compelx variant] -----------------------------------------

template <typename MA, typename VPIV, typename MB>
typename GeMatrix<MA>::IndexType
sv_impl(GeMatrix<MA> &A, DenseVector<VPIV> &piv, GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType IndexType;

    IndexType info = 0;

//
//  Compute the LU factorization of A.
//
    info = trf(A, piv);
    if (info==0) {
//
//      Solve the system A*X = B, overwriting B with X.
//
        trs(NoTrans, A, piv, B);
    }
    return info;
}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (ge)sv [real and complex variant] -----------------------------------------

template <typename MA, typename VPIV, typename MB>
typename GeMatrix<MA>::IndexType
sv_impl(GeMatrix<MA> &A, DenseVector<VPIV> &piv, GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::gesv<IndexType>(A.numRows(),
                                                B.numCols(),
                                                A.data(),
                                                A.leadingDimension(),
                                                piv.data(),
                                                B.data(),
                                                B.leadingDimension());
    ASSERT(info>=0);
    return info;
}

//-- (he)sv [complex variant] -----------------------------------------

template <typename MA, typename VPIV, typename MB, typename VWORK>
typename HeMatrix<MA>::IndexType
sv_impl(HeMatrix<MA> &A, DenseVector<VPIV> &piv, GeMatrix<MB> &B, 
        DenseVector<VWORK> &work)
{
    typedef typename HeMatrix<MA>::IndexType    IndexType;
    typedef typename HeMatrix<MA>::ElementType  ElementType;
    
    if (work.length()==0) {
        ElementType     WORK;
        IndexType       LWORK = -1;
        
        cxxlapack::hesv<IndexType>(getF77Char(A.upLo()),
                                   A.dim(),
                                   B.numCols(),
                                   A.data(),
                                   A.leadingDimension(),
                                   piv.data(),
                                   B.data(),
                                   B.leadingDimension(),
                                   &WORK,
                                   LWORK);
        work.resize(cxxblas::real(WORK));        
    }
    
    IndexType info = cxxlapack::hesv<IndexType>(getF77Char(A.upLo()),
                                                A.dim(),
                                                B.numCols(),
                                                A.data(),
                                                A.leadingDimension(),
                                                piv.data(),
                                                B.data(),
                                                B.leadingDimension(),
                                                work.data(),
                                                work.length());
    ASSERT(info>=0);
    return info;
}

//-- (sy)sv [complex variant] -----------------------------------------

template <typename MA, typename VPIV, typename MB, typename VWORK>
typename SyMatrix<MA>::IndexType
sv_impl(SyMatrix<MA> &A, DenseVector<VPIV> &piv, GeMatrix<MB> &B, 
        DenseVector<VWORK> &work)
{
    typedef typename SyMatrix<MA>::IndexType    IndexType;
    typedef typename SyMatrix<MA>::ElementType  ElementType;
    
    if (work.length()==0) {
        ElementType     WORK;
        IndexType       LWORK = -1;
        
        cxxlapack::sysv<IndexType>(getF77Char(A.upLo()),
                                   A.dim(),
                                   B.numCols(),
                                   A.data(),
                                   A.leadingDimension(),
                                   piv.data(),
                                   B.data(),
                                   B.leadingDimension(),
                                   &WORK,
                                   LWORK);
        work.resize(cxxblas::real(WORK));        
    }


    IndexType info = cxxlapack::sysv<IndexType>(getF77Char(A.upLo()),
                                                A.dim(),
                                                B.numCols(),
                                                A.data(),
                                                A.leadingDimension(),
                                                piv.data(),
                                                B.data(),
                                                B.leadingDimension(),
                                                work.data(),
                                                work.length());
    ASSERT(info>=0);
    return info;
}

//-- (gb)sv [real and complex variant] -----------------------------------------

template <typename MA, typename VPIV, typename MB>
typename GbMatrix<MA>::IndexType
sv_impl(GbMatrix<MA> &A, DenseVector<VPIV> &piv, GeMatrix<MB> &B)
{
    typedef typename GbMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::gbsv<IndexType>(A.numRows(),
                                                A.numSubDiags(),
                                                A.numSuperDiags()-A.numSubDiags(),
                                                B.numCols(),
                                                A.data(),
                                                A.leadingDimension(),
                                                piv.data(),
                                                B.data(),
                                                B.leadingDimension());
    ASSERT(info>=0);
    return info;
}

//-- (sp)sv [complex variant] -----------------------------------------

template <typename MA, typename VPIV, typename MB>
typename HpMatrix<MA>::IndexType
sv_impl(HpMatrix<MA> &A, DenseVector<VPIV> &piv, GeMatrix<MB> &B)
{
    typedef typename HpMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::hpsv<IndexType>(getF77Char(A.upLo()),
                                                A.dim(),
                                                B.numCols(),
                                                A.data(),
                                                piv.data(),
                                                B.data(),
                                                B.leadingDimension());
    ASSERT(info>=0);
    return info;
}

//-- (sp)sv [real variant] -----------------------------------------

template <typename MA, typename VPIV, typename MB>
typename SpMatrix<MA>::IndexType
sv_impl(SpMatrix<MA> &A, DenseVector<VPIV> &piv, GeMatrix<MB> &B)
{
    typedef typename SpMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::spsv<IndexType>(getF77Char(A.upLo()),
                                                A.dim(),
                                                B.numCols(),
                                                A.data(),
                                                piv.data(),
                                                B.data(),
                                                B.leadingDimension());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//-- (ge)sv [real and complex variant] -----------------------------------------

template <typename MA, typename VPIV, typename MB>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsGeMatrix<MB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
sv(MA &&A, VPIV &&piv, MB &&B)
{
    LAPACK_DEBUG_OUT("(ge)sv [real/complex]");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename RemoveRef<VPIV>::Type  VectorPiv;
    typedef typename RemoveRef<MB>::Type    MatrixB;
 
    if (piv.length()<A.numRows()) {
        piv.resize(A.numRows());
    }
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());
    ASSERT((piv.inc()>0 && piv.firstIndex()==1)
        || (piv.inc()<0 && piv.firstIndex()==A.numRows()));
    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==A.numRows());
#   endif
    
#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename MatrixA::NoView    A_org   = A;
    typename VectorPiv::NoView  piv_org = piv;
    typename MatrixB::NoView    B_org   = B;
    
#   endif    
//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::sv_impl(A, piv, B);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename MatrixA::NoView    A_generic   = A;
    typename VectorPiv::NoView  piv_generic = piv;
    typename MatrixB::NoView    B_generic   = B;

    A   = A_org;
    piv = piv_org;
    B   = B_org;

    IndexType _info = external::sv_impl(A, piv, B);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "A_org = " << A_org << std::endl;
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(piv_generic, piv, "piv_generic", "piv")) {
        std::cerr << "piv_org = " << piv_org << std::endl;
        std::cerr << "CXXLAPACK: piv_generic = " << piv_generic << std::endl;
        std::cerr << "F77LAPACK: piv = " << piv << std::endl;
        failed = true;
    }

    if (! isIdentical(B_generic, B, "B_generic", "B")) {
        std::cerr << "B_org = " << B_org << std::endl;
        std::cerr << "CXXLAPACK: B_generic = " << B_generic << std::endl;
        std::cerr << "F77LAPACK: B = " << B << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, " info", "_info")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }
    if (failed) {
        ASSERT(0);
    }

#   endif

    return info;
}

//-- (ge)sv [variant if rhs is vector] -----------------------------------------

template <typename MA, typename VPIV, typename VB>
typename RestrictTo<(IsGeMatrix<MA>::value ||
                     IsHeMatrix<MA>::value ||
                     IsSyMatrix<MA>::value ||
                     IsGbMatrix<MA>::value || 
                     IsHpMatrix<MA>::value || 
                     IsSpMatrix<MA>::value)
                 && IsIntegerDenseVector<VPIV>::value
                 && IsDenseVector<VB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
sv(MA &&A, VPIV &&piv, VB &&b)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VB>::Type    VectorB;

    typedef typename VectorB::ElementType  ElementType;
    typedef typename VectorB::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = MatrixA::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b, n);

    return sv(A, piv, B);
}


#ifdef USE_CXXLAPACK


//-- (he)sv [complex variant] -----------------------------------------
template <typename MA, typename VPIV, typename MB, typename VWORK>
typename RestrictTo<IsHeMatrix<MA>::value
                  && IsIntegerDenseVector<VPIV>::value
                  && IsComplexGeMatrix<MB>::value
                  && IsComplexDenseVector<VWORK>::value,
          typename RemoveRef<MA>::Type::IndexType>::Type
sv(MA &&A, VPIV &&piv, MB &&B, VWORK && work)
{
    LAPACK_DEBUG_OUT("(he)sv [complex]");
    
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    
    if (piv.length()<A.dim()) {
        piv.resize(A.dim());
    }
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT((piv.inc()>0 && piv.firstIndex()==1)
        || (piv.inc()<0 && piv.firstIndex()==A.dim()));
    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==A.dim());
#   endif
//
//  Call implementation
//
    IndexType info = external::sv_impl(A, piv, B, work);

    return info;    
}

//-- (he)sv [complex variant, if rhs is vector] ------------------------
template <typename MA, typename VPIV, typename VB, typename VWORK>
typename RestrictTo<IsHeMatrix<MA>::value
                  && IsIntegerDenseVector<VPIV>::value
                  && IsComplexDenseVector<VB>::value
                  && IsComplexDenseVector<VWORK>::value,
          typename RemoveRef<MA>::Type::IndexType>::Type
sv(MA &&A, VPIV &&piv, VB &&b, VWORK && work)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VB>::Type    VectorB;

    typedef typename VectorB::ElementType  ElementType;
    typedef typename VectorB::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = MatrixA::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b, n);
    
    return sv(A, piv, B, work);
}

//-- (he)sv [complex variant with temporary workspace] ----------------
template <typename MA, typename VPIV, typename MB>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsComplexGeMatrix<MB>::value ,
         typename RemoveRef<MA>::Type::IndexType>::Type
sv(MA &&A, VPIV &&piv, MB &&B)
{
    typedef typename RemoveRef<MA>::Type::Vector WorkVector;

    WorkVector  work;
    
    return sv(A, piv, B, work);
}
    
//-- (sy)sv [real and complex variant] -------------------------------
template <typename MA, typename VPIV, typename MB, typename VWORK>
typename RestrictTo<IsSyMatrix<MA>::value
                  && IsIntegerDenseVector<VPIV>::value
                  && IsGeMatrix<MB>::value
                  && IsDenseVector<VWORK>::value,
          typename RemoveRef<MA>::Type::IndexType>::Type
sv(MA &&A, VPIV &&piv, MB &&B, VWORK && work)
{
    LAPACK_DEBUG_OUT("(sy)sv [real/complex]");
    
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

    if (piv.length()<A.dim()) {
        piv.resize(A.dim());
    }
    
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT((piv.inc()>0 && piv.firstIndex()==1)
        || (piv.inc()<0 && piv.firstIndex()==A.dim()));
    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==A.dim());
#   endif
//
//  Call implementation
//
    IndexType info = external::sv_impl(A, piv, B, work);

    return info;    
}

//-- (sy)sv [real and complex variant, rhs is vector] ---------------------------
template <typename MA, typename VPIV, typename VB, typename VWORK>
typename RestrictTo<IsSyMatrix<MA>::value
                  && IsIntegerDenseVector<VPIV>::value
                  && IsDenseVector<VB>::value
                  && IsDenseVector<VWORK>::value,
          typename RemoveRef<MA>::Type::IndexType>::Type
sv(MA &&A, VPIV &&piv, VB &&b, VWORK && work)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VB>::Type    VectorB;

    typedef typename VectorB::ElementType  ElementType;
    typedef typename VectorB::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = MatrixA::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b, n);
    
    return sv(A, piv, B, work);
}

//-- (sy)sv [real and complex variant with temporary workspace] -----------------

template <typename MA, typename VPIV, typename MB>
typename RestrictTo<IsSyMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsGeMatrix<MB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
sv(MA &&A, VPIV &&piv, MB &&B)
{
    typedef typename RemoveRef<MA>::Type::Vector WorkVector;

    WorkVector  work;
    
    return sv(A, piv, B, work);
}

//-- (gb)sv [real and complex variant] -----------------------------------------

template <typename MA, typename VPIV, typename MB>
typename RestrictTo<IsGbMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsGeMatrix<MB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
sv(MA &&A, VPIV &&piv, MB &&B)
{
    LAPACK_DEBUG_OUT("(gb)sv [real/complex]");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

    if (piv.length()<A.numRows()) {
        piv.resize(A.numRows());
    }
    
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());
    ASSERT((piv.inc()>0 && piv.firstIndex()==1)
        || (piv.inc()<0 && piv.firstIndex()==A.numRows()));
    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==A.numRows());
#   endif
//
//  Call implementation
//
    IndexType info = external::sv_impl(A, piv, B);

    return info;
}

//-- (hp)sv [complex variant] -----------------------------------------

template <typename MA, typename VPIV, typename MB>
typename RestrictTo<IsHpMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsComplexGeMatrix<MB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
sv(MA &&A, VPIV &&piv, MB &&B)
{
    LAPACK_DEBUG_OUT("(hp)sv [complex]");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

    if (piv.length()<A.dim()) {
        piv.resize(A.dim());
    }
    
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstIndex()==1);
    ASSERT((piv.inc()>0 && piv.firstIndex()==1)
        || (piv.inc()<0 && piv.firstIndex()==A.dim()));
    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==A.dim());
#   endif
//
//  Call implementation
//
    IndexType info = external::sv_impl(A, piv, B);

    return info;
}

//-- (sp)sv [real/complex variant] ---------------------------------

template <typename MA, typename VPIV, typename MB>
typename RestrictTo<IsSpMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsGeMatrix<MB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
sv(MA &&A, VPIV &&piv, MB &&B)
{
    LAPACK_DEBUG_OUT("(sp)sv [real/complex]");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

    if (piv.length()<A.dim()) {
        piv.resize(A.dim());
    }
    
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstIndex()==1);
    ASSERT((piv.inc()>0 && piv.firstIndex()==1)
        || (piv.inc()<0 && piv.firstIndex()==A.dim()));
    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==A.dim());
#   endif
//
//  Call implementation
//
    IndexType info = external::sv_impl(A, piv, B);

    return info;
}

#endif // USE_CXXLAPACK
} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_SV_TCC
