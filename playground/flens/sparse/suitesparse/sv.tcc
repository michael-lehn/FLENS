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

#ifndef PLAYGROUND_FLENS_SPARSE_SUITESPARSE_SV_TCC
#define PLAYGROUND_FLENS_SPARSE_SUITESPARSE_SV_TCC 1

#ifdef WITH_UMFPACK

#include <umfpack.h>

namespace flens { namespace suitesparse {

// Interface for AX=B (real)
template <typename MA, typename MX, typename MB>
typename
RestrictTo<IsRealGeCCSMatrix<MA>::value &&
            IsRealGeMatrix<MX>::value &&
            IsRealGeMatrix<MB>::value,
            void>::Type
sv(MA  &&A,
    MX  &&X,
    const MB  &B)
{

    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

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

    const IndexType firstCol = A.firstCol();
    const IndexType firstRow = A.firstRow();

    // SuiteSparse needs IndexBase 0
    // -> shift base temporarily
    if (firstCol != IndexType(0)) {
        A.engine().cols() -= firstCol;
    }
    if (firstRow != IndexType(0)) {
        A.engine().rows() -= firstRow;
    }

    void *Symbolic, *Numeric;

    // symbolic analysis
    umfpack_di_symbolic(A.numRows(), A.numCols(),
                                A.engine().cols().data(),
                                A.engine().rows().data(),
                                A.engine().values().data(),
                                &Symbolic,
                                NULL, NULL);

    // LU factorization
    umfpack_di_numeric(A.engine().cols().data(),
                                A.engine().rows().data(),
                                A.engine().values().data(),
                                Symbolic,
                                &Numeric,
                                NULL, NULL);

    umfpack_di_free_symbolic(&Symbolic);

    // solve system
    for (IndexType i=0; i<B.numCols(); ++i) {
        umfpack_di_solve(UMFPACK_A,
                                  A.engine().cols().data(),
                                  A.engine().rows().data(),
                                  A.engine().values().data(),
                                  X.data() + i*X.leadingDimension(),
                                  B.data() + i*B.leadingDimension(),
                                  Numeric,
                                  NULL, NULL);
    }

    umfpack_di_free_numeric(&Numeric);

    // Reset base to original value
    if (firstCol != IndexType(0)) {
        A.engine().cols() += firstCol;
    }
    if (firstRow != IndexType(0)) {
        A.engine().rows() += firstRow;
    }
}

// Interface for AX=B (complex)
template <typename MA, typename MX, typename MB>
typename
RestrictTo<IsComplexGeCCSMatrix<MA>::value &&
            IsComplexGeMatrix<MX>::value &&
            IsComplexGeMatrix<MB>::value,
            void>::Type
sv(MA  &&A,
    MX  &&X,
    const MB  &B)
{

    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

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

    const IndexType firstCol = A.firstCol();
    const IndexType firstRow = A.firstRow();

    // SuiteSparse needs IndexBase 0
    // -> shift base temporarily
    if (firstCol != IndexType(0)) {
        A.engine().cols() -= firstCol;
    }
    if (firstRow != IndexType(0)) {
        A.engine().rows() -= firstRow;
    }

    void *Symbolic, *Numeric;

    // symbolic analysis
    umfpack_zi_symbolic(A.numRows(), A.numCols(),
                                A.engine().cols().data(),
                                A.engine().rows().data(),
                                reinterpret_cast<double*>(A.engine().values().data()), NULL,
                                &Symbolic,
                                NULL, NULL);

    // LU factorization
    umfpack_zi_numeric(A.engine().cols().data(),
                                A.engine().rows().data(),
                                reinterpret_cast<double*>(A.engine().values().data()), NULL,
                                Symbolic,
                                &Numeric,
                                NULL, NULL);
    umfpack_zi_free_symbolic(&Symbolic);

    for (IndexType i=0; i<B.numCols(); ++i) {
        // solve system
        umfpack_zi_solve(UMFPACK_A,
                                  A.engine().cols().data(),
                                  A.engine().rows().data(),
                                  reinterpret_cast<double*>(A.engine().values().data()), NULL,
                                  reinterpret_cast<double*>(X.data()+i*X.leadingDimension()), NULL,
                                  reinterpret_cast<const double*>(B.data()+i*B.leadingDimension()), NULL,
                                  Numeric,
                                  NULL, NULL);
    }
    umfpack_zi_free_numeric(&Numeric);

    // Reset base to original value
    if (firstCol != IndexType(0)) {
        A.engine().cols() += firstCol;
    }
    if (firstRow != IndexType(0)) {
        A.engine().rows() += firstRow;
    }
}

// Interface for vectors
template <typename MA, typename VX, typename VB>
typename
RestrictTo<IsGeCCSMatrix<MA>::value &&
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

} } // namespace suitesparse, flens

#endif // WITH_UMFPACK

#endif // PLAYGROUND_FLENS_SPARSE_SUITESPARSE_SV_TCC
