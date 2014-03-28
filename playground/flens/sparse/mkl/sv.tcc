/*
 *   Copyright (c) 2013, Klaus Pototzky
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

#ifndef PLAYGROUND_FLENS_SPARSE_MKL_SV_TCC
#define PLAYGROUND_FLENS_SPARSE_MKL_SV_TCC 1

#ifdef WITH_MKLDSS

#include <mkl_dss.h>
#include <mkl_types.h>

#include<flens/auxiliary/auxiliary.h>
#include<flens/matrixtypes/matrixtypes.h>
#include<flens/vectortypes/vectortypes.h>

namespace flens { namespace mkldss {

template <typename MA, typename MX, typename MB>
typename
RestrictTo<IsRealGeCCSMatrix<MA>::value &&
           IsRealGeMatrix<MX>::value &&
           IsRealGeMatrix<MB>::value,
           void>::Type
sv(Transpose transA,
   const MA  &A,
   MX        &&X,
   const MB  &B)
{
    transA = Transpose(transA^Trans);
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename MatrixA::ElementType   ElementType;

    ASSERT(X.order()==ColMajor);
    ASSERT(B.order()==ColMajor);

    if (X.numCols()==IndexType(0) || X.numRows()==IndexType(0)) {
        X.resize(B.numRows(), B.numRows(),
                 B.firstRow(), B.firstCol());
    }

    ASSERT(X.numRows()==B.numRows());
    ASSERT(X.numCols()==B.numCols());
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.numCols()==B.numRows());

    const IndexType firstIndex = A.firstCol();
    ASSERT(A.firstRow()==A.firstCol());
    ASSERT(firstIndex==IndexType(0) ||
           firstIndex==IndexType(1));


    /* Allocate storage for the solver handle and the right-hand side. */
    void *handle;
    int  error;
    IndexType opt  = MKL_DSS_DEFAULTS;
    IndexType sym  = MKL_DSS_NON_SYMMETRIC;
    IndexType type = MKL_DSS_INDEFINITE;

    IndexType nRows     = A.numRows();
    IndexType nCols     = A.numCols();
    IndexType nNonZeros = A.engine().values().length();

    // Set flag for single precision
    if (IsSame<ElementType, float>::value) {
        opt |= MKL_DSS_SINGLE_PRECISION;
    }

    // Set flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt |= MKL_DSS_ZERO_BASED_INDEXING;
    }

    error = dss_create (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);

    // Unset flag for single precision
    if (IsSame<ElementType, float>::value) {
        opt &= ~MKL_DSS_SINGLE_PRECISION;
    }
    // Unset flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt &= ~MKL_DSS_ZERO_BASED_INDEXING;
    }

    error = dss_define_structure (handle,
                                  sym,
                                  A.engine().cols().data(),
                                  nRows, nCols,
                                  A.engine().rows().data(), nNonZeros);

    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_reorder (handle, opt, 0);
    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_factor_real (handle, type, A.engine().values().data());
    ASSERT(error==MKL_DSS_SUCCESS);

    // Add Transpose flag
    if (transA==Trans || transA==ConjTrans) {
        opt |= MKL_DSS_TRANSPOSE_SOLVE;

    }

    IndexType nRhs = B.numCols();
    if (B.leadingDimension()!=B.numRows()
     || X.leadingDimension() != X.numRows())
    {
        error = dss_solve_real (handle, opt, B.data(), nRhs, X.data());
        ASSERT(error==MKL_DSS_SUCCESS);
    } else {
        IndexType one(1);
        for (IndexType i=0; i<nRhs; ++i) {
            error = dss_solve_real(handle, opt, B.data()+i*B.leadingDimension(),
                                   one, X.data()+i*X.leadingDimension());
            ASSERT(error==MKL_DSS_SUCCESS);
        }
    }

    // Remove Transpose flag
    if (transA==Trans || transA==ConjTrans) {
        opt &= ~MKL_DSS_TRANSPOSE_SOLVE;

    }

    error = dss_delete (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);
}

template <typename MA, typename MX, typename MB>
typename
RestrictTo<IsRealGeCRSMatrix<MA>::value &&
           IsRealGeMatrix<MX>::value &&
           IsRealGeMatrix<MB>::value,
           void>::Type
sv(Transpose transA,
   const MA        &A,
   MX        &&X,
   const MB  &B)
{
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename MatrixA::ElementType   ElementType;

    ASSERT(X.order()==ColMajor);
    ASSERT(B.order()==ColMajor);

    if (X.numCols()==IndexType(0) || X.numRows()==IndexType(0)) {
        X.resize(B.numRows(), B.numRows(),
                 B.firstRow(), B.firstCol());
    }

    ASSERT(X.numRows()==B.numRows());
    ASSERT(X.numCols()==B.numCols());
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.numCols()==B.numRows());

    const IndexType firstIndex = A.firstCol();
    ASSERT(A.firstRow()==A.firstCol());
    ASSERT(firstIndex==IndexType(0) ||
           firstIndex==IndexType(1));


    /* Allocate storage for the solver handle and the right-hand side. */
    void *handle;
    int  error;
    IndexType opt  = MKL_DSS_DEFAULTS;
    IndexType sym  = MKL_DSS_NON_SYMMETRIC;
    IndexType type = MKL_DSS_INDEFINITE;

    IndexType nRows     = A.numRows();
    IndexType nCols     = A.numCols();
    IndexType nNonZeros = A.engine().values().length();

    // Set flag for single precision
    if (IsSame<ElementType, float>::value) {
        opt |= MKL_DSS_SINGLE_PRECISION;
    }

    // Set flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt |= MKL_DSS_ZERO_BASED_INDEXING;
    }

    error = dss_create (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);

    // Unset flag for single precision
    if (IsSame<ElementType, float>::value) {
        opt &= ~MKL_DSS_SINGLE_PRECISION;
    }
    // Unset flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt &= ~MKL_DSS_ZERO_BASED_INDEXING;
    }

    error = dss_define_structure (handle,
                                  sym,
                                  A.engine().rows().data(),
                                  nRows, nCols,
                                  A.engine().cols().data(), nNonZeros);

    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_reorder (handle, opt, 0);
    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_factor_real (handle, type, A.engine().values().data());
    ASSERT(error==MKL_DSS_SUCCESS);

    // Add Transpose flag
    if (transA==Trans || transA==ConjTrans) {
        opt |= MKL_DSS_TRANSPOSE_SOLVE;

    }

    IndexType nRhs = B.numCols();
    if (B.leadingDimension()!=B.numRows() || X.leadingDimension() != X.numRows()) {

        error = dss_solve_real (handle, opt, B.data(), nRhs, X.data());
        ASSERT(error==MKL_DSS_SUCCESS);

    } else {

        IndexType one(1);
        for (IndexType i=0; i<nRhs; ++i) {
            error = dss_solve_real (handle, opt, B.data()+i*B.leadingDimension(), one, X.data()+i*X.leadingDimension());
            ASSERT(error==MKL_DSS_SUCCESS);
        }
    }

    // Remove Transpose flag
    if (transA==Trans || transA==ConjTrans) {
        opt &= ~MKL_DSS_TRANSPOSE_SOLVE;

    }

    error = dss_delete (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);

}

template <typename MA, typename MX, typename MB>
typename
RestrictTo<IsComplexGeCCSMatrix<MA>::value &&
           IsComplexGeMatrix<MX>::value &&
           IsComplexGeMatrix<MB>::value,
           void>::Type
sv(Transpose transA,
   MA        &&A,
   MX        &&X,
   const MB  &B)
{
    transA = Transpose(transA^Trans); 
    
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename MatrixA::ElementType     ElementType;

    ASSERT(X.order()==ColMajor);
    ASSERT(B.order()==ColMajor);
    ASSERT(transA!=Conj);

    if (X.numCols()==IndexType(0) || X.numRows()==IndexType(0)) {
        X.resize(B.numRows(), B.numRows(),
                 B.firstRow(), B.firstCol());
    }

    ASSERT(X.numRows()==B.numRows());
    ASSERT(X.numCols()==B.numCols());
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.numCols()==B.numRows());

    const IndexType firstIndex = A.firstCol();
    ASSERT(A.firstRow()==A.firstCol());
    ASSERT(firstIndex==IndexType(0) ||
           firstIndex==IndexType(1));


    /* Allocate storage for the solver handle and the right-hand side. */
    void *handle;
    int  error;
    IndexType opt  = MKL_DSS_DEFAULTS;
    IndexType sym  = MKL_DSS_NON_SYMMETRIC;
    IndexType type = MKL_DSS_INDEFINITE;

    IndexType nRows     = A.numRows();
    IndexType nCols     = A.numCols();
    IndexType nNonZeros = A.engine().values().length();


    // Set flag for single precision
    if (IsSame<ElementType, std::complex<float> >::value) {
        opt |= MKL_DSS_SINGLE_PRECISION;
    }

    // Set flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt |= MKL_DSS_ZERO_BASED_INDEXING;
    }

    error = dss_create (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);

    // Unset flag for single precision
    if (IsSame<ElementType, std::complex<float> >::value) {
        opt &= ~MKL_DSS_SINGLE_PRECISION;
    }
    // Unset flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt &= ~MKL_DSS_ZERO_BASED_INDEXING;
    }

    error = dss_define_structure (handle,
                                  sym,
                                  A.engine().cols().data(),
                                  nRows, nCols,
                                  A.engine().rows().data(), nNonZeros);

    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_reorder (handle, opt, 0);
    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_factor_complex (handle, type, A.engine().values().data());
    ASSERT(error==MKL_DSS_SUCCESS);

    // Add Transpose flag
    if (transA==Trans) {
        opt |= MKL_DSS_TRANSPOSE_SOLVE;
    } else if (transA==ConjTrans) {
        opt |= MKL_DSS_CONJUGATE_SOLVE;

    }

    IndexType nRhs = B.numCols();
    if (B.leadingDimension()!=B.numRows() || X.leadingDimension() != X.numRows()) {

        error = dss_solve_complex (handle, opt, B.data(), nRhs, X.data());
        ASSERT(error==MKL_DSS_SUCCESS);

    } else {

        IndexType one(1);
        for (IndexType i=0; i<nRhs; ++i) {
            error = dss_solve_complex (handle, opt, B.data()+i*B.leadingDimension(), one, X.data()+i*X.leadingDimension());
            ASSERT(error==MKL_DSS_SUCCESS);
        }
    }

    // Remove Transpose flag
    if (transA==Trans) {
        opt &= ~MKL_DSS_TRANSPOSE_SOLVE;
    } else if (transA==ConjTrans) {
        opt &= ~MKL_DSS_CONJUGATE_SOLVE;
    }

    error = dss_delete (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);
}

template <typename MA, typename MX, typename MB>
typename
RestrictTo<IsComplexGeCRSMatrix<MA>::value &&
           IsComplexGeMatrix<MX>::value &&
           IsComplexGeMatrix<MB>::value,
           void>::Type
sv(Transpose transA,
   MA        &&A,
   MX        &&X,
   const MB  &B)
{
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename MatrixA::ElementType     ElementType;

    ASSERT(X.order()==ColMajor);
    ASSERT(B.order()==ColMajor);
    ASSERT(transA!=Conj);

    if (X.numCols()==IndexType(0) || X.numRows()==IndexType(0)) {
        X.resize(B.numRows(), B.numRows(),
                 B.firstRow(), B.firstCol());
    }

    ASSERT(X.numRows()==B.numRows());
    ASSERT(X.numCols()==B.numCols());
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.numCols()==B.numRows());

    const IndexType firstIndex = A.firstCol();
    ASSERT(A.firstRow()==A.firstCol());
    ASSERT(firstIndex==IndexType(0) ||
           firstIndex==IndexType(1));


    /* Allocate storage for the solver handle and the right-hand side. */
    void *handle;
    int  error;
    IndexType opt  = MKL_DSS_DEFAULTS;
    IndexType sym  = MKL_DSS_NON_SYMMETRIC;
    IndexType type = MKL_DSS_INDEFINITE;

    IndexType nRows     = A.numRows();
    IndexType nCols     = A.numCols();
    IndexType nNonZeros = A.engine().values().length();


    // Set flag for single precision
    if (IsSame<ElementType, std::complex<float> >::value) {
        opt |= MKL_DSS_SINGLE_PRECISION;
    }

    // Set flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt |= MKL_DSS_ZERO_BASED_INDEXING;
    }

    error = dss_create (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);

    // Unset flag for single precision
    if (IsSame<ElementType, std::complex<float> >::value) {
        opt &= ~MKL_DSS_SINGLE_PRECISION;
    }
    // Unset flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt &= ~MKL_DSS_ZERO_BASED_INDEXING;
    }

    error = dss_define_structure (handle,
                                  sym,
                                  A.engine().rows().data(),
                                  nRows, nCols,
                                  A.engine().cols().data(), nNonZeros);

    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_reorder (handle, opt, 0);
    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_factor_complex (handle, type, A.engine().values().data());
    ASSERT(error==MKL_DSS_SUCCESS);

    // Add Transpose flag
    if (transA==Trans) {
        opt |= MKL_DSS_TRANSPOSE_SOLVE;
    } else if (transA==ConjTrans) {
        opt |= MKL_DSS_CONJUGATE_SOLVE;

    }

    IndexType nRhs = B.numCols();
    if (B.leadingDimension()!=B.numRows() || X.leadingDimension() != X.numRows()) {

        error = dss_solve_complex (handle, opt, B.data(), nRhs, X.data());
        ASSERT(error==MKL_DSS_SUCCESS);

    } else {

        IndexType one(1);
        for (IndexType i=0; i<nRhs; ++i) {
            error = dss_solve_complex (handle, opt, B.data()+i*B.leadingDimension(), one, X.data()+i*X.leadingDimension());
            ASSERT(error==MKL_DSS_SUCCESS);
        }
    }

    // Remove Transpose flag
    if (transA==Trans) {
        opt &= ~MKL_DSS_TRANSPOSE_SOLVE;
    } else if (transA==ConjTrans) {
        opt &= ~MKL_DSS_CONJUGATE_SOLVE;
    }

    error = dss_delete (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);
}


template <typename MA, typename MX, typename MB>
typename
RestrictTo<IsRealSyCCSMatrix<MA>::value &&
           IsRealGeMatrix<MX>::value &&
           IsRealGeMatrix<MB>::value,
           void>::Type
sv(MA        &&A,
   MX        &&X,
   const MB  &B)
{
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename MatrixA::ElementType   ElementType;

    ASSERT(A.upLo()==Lower);
    ASSERT(X.order()==ColMajor);
    ASSERT(B.order()==ColMajor);

    if (X.numCols()==IndexType(0) || X.numRows()==IndexType(0)) {
        X.resize(B.numRows(), B.numRows(),
                 B.firstRow(), B.firstCol());
    }

    ASSERT(X.numRows()==B.numRows());
    ASSERT(X.numCols()==B.numCols());
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.numCols()==B.numRows());

    const IndexType firstIndex = A.firstCol();
    ASSERT(A.firstRow()==A.firstCol());
    ASSERT(firstIndex==IndexType(0) ||
           firstIndex==IndexType(1));


    /* Allocate storage for the solver handle and the right-hand side. */
    void *handle;
    int  error;
    IndexType opt  = MKL_DSS_DEFAULTS;
    IndexType sym  = MKL_DSS_SYMMETRIC;
    IndexType type = MKL_DSS_INDEFINITE;

    IndexType nRows     = A.numRows();
    IndexType nCols     = A.numCols();
    IndexType nNonZeros = A.engine().values().length();

    // Set flag for single precision
    if (IsSame<ElementType, float>::value) {
        opt |= MKL_DSS_SINGLE_PRECISION;
    }

    // Set flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt |= MKL_DSS_ZERO_BASED_INDEXING;
    }

    error = dss_create (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);

    // Unset flag for single precision
    if (IsSame<ElementType, float>::value) {
        opt &= ~MKL_DSS_SINGLE_PRECISION;
    }
    // Unset flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt &= ~MKL_DSS_ZERO_BASED_INDEXING;
    }

    error = dss_define_structure (handle,
                                  sym,
                                  A.engine().cols().data(),
                                  nRows, nCols,
                                  A.engine().rows().data(), nNonZeros);

    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_reorder (handle, opt, 0);
    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_factor_real (handle, type, A.engine().values().data());
    ASSERT(error==MKL_DSS_SUCCESS);

    IndexType nRhs = B.numCols();
    if (B.leadingDimension()!=B.numRows() || X.leadingDimension() != X.numRows()) {

        error = dss_solve_real (handle, opt, B.data(), nRhs, X.data());
        ASSERT(error==MKL_DSS_SUCCESS);

    } else {

        IndexType one(1);
        for (IndexType i=0; i<nRhs; ++i) {
            error = dss_solve_real (handle, opt, B.data()+i*B.leadingDimension(), one, X.data()+i*X.leadingDimension());
            ASSERT(error==MKL_DSS_SUCCESS);
        }
    }

    error = dss_delete (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);

}


template <typename MA, typename MX, typename MB>
typename
RestrictTo<IsRealSyCRSMatrix<MA>::value &&
           IsRealGeMatrix<MX>::value &&
           IsRealGeMatrix<MB>::value,
           void>::Type
sv(MA        &&A,
   MX        &&X,
   const MB  &B)
{
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename MatrixA::ElementType   ElementType;

    ASSERT(A.upLo()==Upper);
    ASSERT(X.order()==ColMajor);
    ASSERT(B.order()==ColMajor);

    if (X.numCols()==IndexType(0) || X.numRows()==IndexType(0)) {
        X.resize(B.numRows(), B.numRows(),
                 B.firstRow(), B.firstCol());
    }

    ASSERT(X.numRows()==B.numRows());
    ASSERT(X.numCols()==B.numCols());
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.numCols()==B.numRows());

    const IndexType firstIndex = A.firstCol();
    ASSERT(A.firstRow()==A.firstCol());
    ASSERT(firstIndex==IndexType(0) ||
           firstIndex==IndexType(1));


    /* Allocate storage for the solver handle and the right-hand side. */
    void *handle;
    int  error;
    IndexType opt  = MKL_DSS_DEFAULTS;
    IndexType sym  = MKL_DSS_SYMMETRIC;
    IndexType type = MKL_DSS_INDEFINITE;

    IndexType nRows     = A.numRows();
    IndexType nCols     = A.numCols();
    IndexType nNonZeros = A.engine().values().length();

    // Set flag for single precision
    if (IsSame<ElementType, float>::value) {
        opt |= MKL_DSS_SINGLE_PRECISION;
    }

    // Set flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt |= MKL_DSS_ZERO_BASED_INDEXING;
    }

    error = dss_create (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);

    // Unset flag for single precision
    if (IsSame<ElementType, float>::value) {
        opt &= ~MKL_DSS_SINGLE_PRECISION;
    }
    // Unset flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt &= ~MKL_DSS_ZERO_BASED_INDEXING;
    }

    error = dss_define_structure (handle,
                                  sym,
                                  A.engine().rows().data(),
                                  nRows, nCols,
                                  A.engine().cols().data(), nNonZeros);

    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_reorder (handle, opt, 0);
    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_factor_real (handle, type, A.engine().values().data());
    ASSERT(error==MKL_DSS_SUCCESS);

    IndexType nRhs = B.numCols();
    if (B.leadingDimension()!=B.numRows() || X.leadingDimension() != X.numRows()) {

        error = dss_solve_real (handle, opt, B.data(), nRhs, X.data());
        ASSERT(error==MKL_DSS_SUCCESS);

    } else {

        IndexType one(1);
        for (IndexType i=0; i<nRhs; ++i) {
            error = dss_solve_real (handle, opt, B.data()+i*B.leadingDimension(), one, X.data()+i*X.leadingDimension());
            ASSERT(error==MKL_DSS_SUCCESS);
        }
    }

    error = dss_delete (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);

}

template <typename MA, typename MX, typename MB>
typename
RestrictTo<IsComplexSyCCSMatrix<MA>::value &&
           IsComplexGeMatrix<MX>::value &&
           IsComplexGeMatrix<MB>::value,
           void>::Type
sv(MA        &&A,
   MX        &&X,
   const MB  &B)
{
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename MatrixA::ElementType   ElementType;

    ASSERT(A.upLo()==Lower);
    ASSERT(X.order()==ColMajor);
    ASSERT(B.order()==ColMajor);

    if (X.numCols()==IndexType(0) || X.numRows()==IndexType(0)) {
        X.resize(B.numRows(), B.numRows(),
                 B.firstRow(), B.firstCol());
    }

    ASSERT(X.numRows()==B.numRows());
    ASSERT(X.numCols()==B.numCols());
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.numCols()==B.numRows());

    const IndexType firstIndex = A.firstCol();
    ASSERT(A.firstRow()==A.firstCol());
    ASSERT(firstIndex==IndexType(0) ||
           firstIndex==IndexType(1));


    /* Allocate storage for the solver handle and the right-hand side. */
    void *handle;
    int  error;
    IndexType opt  = MKL_DSS_DEFAULTS;
    IndexType sym  = MKL_DSS_SYMMETRIC;
    IndexType type = MKL_DSS_INDEFINITE;

    IndexType nRows     = A.numRows();
    IndexType nCols     = A.numCols();
    IndexType nNonZeros = A.engine().values().length();

    // Set flag for single precision
    if (IsSame<ElementType, std::complex<float> >::value) {
        opt |= MKL_DSS_SINGLE_PRECISION;
    }

    // Set flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt |= MKL_DSS_ZERO_BASED_INDEXING;
    }

    error = dss_create (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);

    // Unset flag for single precision
    if (IsSame<ElementType, std::complex<float> >::value) {
        opt &= ~MKL_DSS_SINGLE_PRECISION;
    }
    // Unset flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt &= ~MKL_DSS_ZERO_BASED_INDEXING;
    }
    error = dss_define_structure (handle,
                                  sym,
                                  A.engine().cols().data(),
                                  nRows, nCols,
                                  A.engine().rows().data(), nNonZeros);

    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_reorder (handle, opt, 0);
    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_factor_complex (handle, type, A.engine().values().data());
    ASSERT(error==MKL_DSS_SUCCESS);

    IndexType nRhs = B.numCols();
    if (B.leadingDimension()!=B.numRows() || X.leadingDimension() != X.numRows()) {

        error = dss_solve_complex (handle, opt, B.data(), nRhs, X.data());
        ASSERT(error==MKL_DSS_SUCCESS);

    } else {

        IndexType one(1);
        for (IndexType i=0; i<nRhs; ++i) {
            error = dss_solve_complex (handle, opt, B.data()+i*B.leadingDimension(), one, X.data()+i*X.leadingDimension());
            ASSERT(error==MKL_DSS_SUCCESS);
        }
    }

    error = dss_delete (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);
}

template <typename MA, typename MX, typename MB>
typename
RestrictTo<IsComplexSyCRSMatrix<MA>::value &&
           IsComplexGeMatrix<MX>::value &&
           IsComplexGeMatrix<MB>::value,
           void>::Type
sv(MA        &&A,
   MX        &&X,
   const MB  &B)
{
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename MatrixA::ElementType   ElementType;

    ASSERT(A.upLo()==Upper);
    ASSERT(X.order()==ColMajor);
    ASSERT(B.order()==ColMajor);

    if (X.numCols()==IndexType(0) || X.numRows()==IndexType(0)) {
        X.resize(B.numRows(), B.numRows(),
                 B.firstRow(), B.firstCol());
    }

    ASSERT(X.numRows()==B.numRows());
    ASSERT(X.numCols()==B.numCols());
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.numCols()==B.numRows());

    const IndexType firstIndex = A.firstCol();
    ASSERT(A.firstRow()==A.firstCol());
    ASSERT(firstIndex==IndexType(0) ||
           firstIndex==IndexType(1));


    /* Allocate storage for the solver handle and the right-hand side. */
    void *handle;
    int  error;
    IndexType opt  = MKL_DSS_DEFAULTS;
    IndexType sym  = MKL_DSS_SYMMETRIC;
    IndexType type = MKL_DSS_INDEFINITE;

    IndexType nRows     = A.numRows();
    IndexType nCols     = A.numCols();
    IndexType nNonZeros = A.engine().values().length();

    // Set flag for single precision
    if (IsSame<ElementType, std::complex<float> >::value) {
        opt |= MKL_DSS_SINGLE_PRECISION;
    }

    // Set flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt |= MKL_DSS_ZERO_BASED_INDEXING;
    }

    error = dss_create (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);

    // Unset flag for single precision
    if (IsSame<ElementType, std::complex<float> >::value) {
        opt &= ~MKL_DSS_SINGLE_PRECISION;
    }
    // Unset flag for IndexBase
    if (firstIndex==IndexType(0)) {
        opt &= ~MKL_DSS_ZERO_BASED_INDEXING;
    }
    error = dss_define_structure (handle,
                                  sym,
                                  A.engine().rows().data(),
                                  nRows, nCols,
                                  A.engine().cols().data(), nNonZeros);

    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_reorder (handle, opt, 0);
    ASSERT(error==MKL_DSS_SUCCESS);

    error = dss_factor_complex (handle, type, A.engine().values().data());
    ASSERT(error==MKL_DSS_SUCCESS);

    IndexType nRhs = B.numCols();
    if (B.leadingDimension()!=B.numRows() || X.leadingDimension() != X.numRows()) {

        error = dss_solve_complex (handle, opt, B.data(), nRhs, X.data());
        ASSERT(error==MKL_DSS_SUCCESS);

    } else {

        IndexType one(1);
        for (IndexType i=0; i<nRhs; ++i) {
            error = dss_solve_complex (handle, opt, B.data()+i*B.leadingDimension(), one, X.data()+i*X.leadingDimension());
            ASSERT(error==MKL_DSS_SUCCESS);
        }
    }

    error = dss_delete (handle, opt);
    ASSERT(error==MKL_DSS_SUCCESS);
}


// Interface for vectors
template <typename MA, typename VX, typename VB>
typename
RestrictTo<(IsGeCCSMatrix<MA>::value || IsGeCRSMatrix<MA>::value) &&
           IsDenseVector<VX>::value &&
           IsDenseVector<VB>::value,
           void>::Type
sv(Transpose transA,
   MA  &&A,
   VX  &&x,
   const VB  &b)
{

    typedef typename RemoveRef<VB>::Type    VectorB;
    typedef typename VectorB::ElementType   ElementType;
    typedef typename VectorB::IndexType     IndexType;

    if (x.length()==IndexType(0)) {
        x.resize(b.length());
    }
    ASSERT(x.length()==b.length());
    ASSERT(b.stride()==IndexType(1));
    ASSERT(x.stride()==IndexType(1));

    IndexType n      = b.length();
    const GeMatrix<ConstFullStorageView<ElementType, ColMajor> >  B(n, 1, b, n);
    GeMatrix<FullStorageView<ElementType, ColMajor> >  X(n, 1, x, n);

    return sv(transA, A, X, B);
}

template <typename MA, typename VX, typename VB>
typename
RestrictTo<(IsSyCCSMatrix<MA>::value || IsSyCRSMatrix<MA>::value) &&
           IsDenseVector<VX>::value &&
           IsDenseVector<VB>::value,
           void>::Type
sv(MA  &&A,
   VX  &&x,
   const VB  &b)
{

    typedef typename RemoveRef<VB>::Type    VectorB;
    typedef typename VectorB::ElementType   ElementType;
    typedef typename VectorB::IndexType     IndexType;

    if (x.length()==IndexType(0)) {
        x.resize(b.length());
    }
    ASSERT(x.length()==b.length());
    ASSERT(b.stride()==IndexType(1));
    ASSERT(x.stride()==IndexType(1));

    IndexType n      = b.length();
    const GeMatrix<ConstFullStorageView<ElementType, ColMajor> >  B(n, 1, b, n);
    GeMatrix<FullStorageView<ElementType, ColMajor> >  X(n, 1, x, n);

    return sv(A, X, B);
}

} } // namespace mkldss, flens

#endif // WITH_MKLDSS

#endif // PLAYGROUND_FLENS_SPARSE_MKL_SV_TCC
