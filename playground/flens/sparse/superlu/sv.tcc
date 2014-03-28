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

#ifndef PLAYGROUND_FLENS_SPARSE_SUPERLU_SV_TCC
#define PLAYGROUND_FLENS_SPARSE_SUPERLU_SV_TCC 1

#ifdef WITH_SUPERLU

namespace flens { namespace superlu {

template <typename MA, typename PC, typename PR, typename MB>
typename
RestrictTo<IsGeCCSMatrix<MA>::value  &&
           IsDenseVector<PC>::value &&
           IsDenseVector<PR>::value &&
           IsGeMatrix<MB>::value,
           int>::Type
sv(MA              &&A,
    PC              &&pc,
    PR              &&pr,
    MB              &&B)
{
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::ElementType   ElementType;
    typedef typename MatrixA::IndexType     IndexType;

    if (pc.length()==0) {
        pc.resize(A.numCols());
    }
    if (pr.length()==0) {
        pr.resize(A.numRows());
    }

    ASSERT(A.numRows()==A.numCols());
    ASSERT(pr.length()==A.numRows());
    ASSERT(pc.length()==A.numCols());

    ASSERT(pc.stride()==1);
    ASSERT(pr.stride()==1);
    ASSERT(B.order()==ColMajor);

    const IndexType firstCol = A.firstCol();
    const IndexType firstRow = A.firstRow();

    // SuperLu needs IndexBase 0
    // -> shift base temporarily
    if (firstCol != IndexType(0)) {
        A.engine().cols() -= firstCol;
    }
    if (firstRow != IndexType(0)) {
        A.engine().rows() -= firstRow;
    }

    superlu_options_t   options;
    SuperLUStat_t       stat;
    SuperMatrix         _A, _L, _U, _B;
    Create_CompCol_Matrix(&_A,
                          A.numRows(), A.numCols(), A.engine().numNonZeros(),
                          A.engine().values().data(),
                          A.engine().rows().data(),
                          A.engine().cols().data(),
                          SLU_NC, SLU_GE);

    Create_Dense_Matrix(&_B, B.numRows(), B.numCols(),
                        B.data(), B.leadingDimension(),
                        SLU_DN, SLU_GE);

    set_default_options(&options);
    options.ColPerm = NATURAL;
    StatInit(&stat);
    int info = gssv<ElementType>(&options, &_A, pc.data(), pr.data(),
                                 &_L, &_U, &_B, &stat);
    Destroy_SuperMatrix_Store(&_A);
    Destroy_SuperMatrix_Store(&_B);
    Destroy_SuperNode_Matrix(&_L);
    Destroy_CompCol_Matrix(&_U);
    StatFree(&stat);

    // Reset base to original value
    if (firstCol != IndexType(0)) {
        A.engine().cols() += firstCol;
    }
    if (firstRow != IndexType(0)) {
        A.engine().rows() += firstRow;
    }

    return info;
}

template <typename MA, typename PC, typename PR, typename MB>
typename
RestrictTo<IsGeCRSMatrix<MA>::value  &&
           IsDenseVector<PC>::value &&
           IsDenseVector<PR>::value &&
           IsGeMatrix<MB>::value,
           int>::Type
sv(MA              &&A,
    PC              &&pc,
    PR              &&pr,
    MB              &&B)
{
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::ElementType   ElementType;
    typedef typename MatrixA::IndexType     IndexType;

    if (pc.length()==0) {
        pc.resize(A.numCols());
    }
    if (pr.length()==0) {
        pr.resize(A.numRows());
    }

    // Square Matrix
    ASSERT(A.numRows()==A.numCols());
    ASSERT(pr.length()==A.numRows());
    ASSERT(pc.length()==A.numCols());

    ASSERT(pc.stride()==1);
    ASSERT(pr.stride()==1);
    ASSERT(B.order()==ColMajor);

    const IndexType firstCol = A.firstCol();
    const IndexType firstRow = A.firstRow();

    // SuperLu needs IndexBase 0
    // -> shift base temporarily

    if (firstCol != IndexType(0)) {
        A.engine().cols() -= firstCol;
    }
    if (firstRow != IndexType(0)) {
        A.engine().rows() -= firstRow;
    }

    superlu_options_t   options;
    SuperLUStat_t       stat;
    SuperMatrix         _A, _L, _U, _B;


    Create_CompCol_Matrix(&_A,
                          A.numCols(), A.numRows(), A.engine().numNonZeros(),
                          A.engine().values().data(),
                          A.engine().cols().data(),
                          A.engine().rows().data(),
                          SLU_NR, SLU_GE);



    Create_Dense_Matrix(&_B, B.numRows(), B.numCols(),
                        B.data(), B.leadingDimension(),
                        SLU_DN, SLU_GE);

    set_default_options(&options);
    options.ColPerm = NATURAL;
    StatInit(&stat);
    int info = gssv<ElementType>(&options, &_A, pc.data(), pr.data(),
                                 &_L, &_U, &_B, &stat);
    Destroy_SuperMatrix_Store(&_A);
    Destroy_SuperMatrix_Store(&_B);
    Destroy_SuperNode_Matrix(&_L);
    Destroy_CompCol_Matrix(&_U);
    StatFree(&stat);

    // Reset base to original value
    if (firstCol != IndexType(0)) {
        A.engine().cols() += firstCol;
    }
    if (firstRow != IndexType(0)) {
        A.engine().rows() += firstRow;
    }

    return info;
}

template <typename MA, typename PC, typename PR, typename VB>
typename
RestrictTo<(IsGeCCSMatrix<MA>::value || IsGeCRSMatrix<MA>::value) &&
            IsDenseVector<PC>::value &&
            IsDenseVector<PR>::value &&
            IsDenseVector<VB>::value,
            int>::Type
sv(MA              &&A,
    PC              &&pc,
    PR              &&pr,
    VB              &&b)
{

    typedef typename RemoveRef<VB>::Type    VectorB;
    typedef typename VectorB::ElementType   ElementType;
    typedef typename VectorB::IndexType     IndexType;

    ASSERT(b.stride()==IndexType(1));
    IndexType n      = b.length();
    GeMatrix<FullStorageView<ElementType, ColMajor> >  B(n, 1, b, n);
    return superlu::sv(A, pc, pr, B);
}

} } // namespace superlu, flens

#endif // WITH_SUPERLU

#endif // PLAYGROUND_FLENS_SPARSE_SUPERLU_SV_TCC
