/*
 *   Copyright (c) 2010, Michael Lehn
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

#include <flens/blas/blaslogclear.h>

#define FLENS_BLASLOG_SETTAG(tag)
#define FLENS_BLASLOG_UNSETTAG

#define FLENS_BLASLOG_END
#define FLENS_BLASLOG_BEGIN_ASSIGNMENT(X, Y)
#define FLENS_BLASLOG_BEGIN_PLUSASSIGNMENT(X, Y)
#define FLENS_BLASLOG_BEGIN_MINUSASSIGNMENT(X, Y)

#define FLENS_BLASLOG_IDENTICAL(A, B)
#define FLENS_BLASLOG_TMP_ADD(tmp)
#define FLENS_BLASLOG_TMP_TRON
#define FLENS_BLASLOG_TMP_TROFF
#define FLENS_BLASLOG_TMP_REMOVE(tmp, X)

#define FLENS_BLASLOG_ERROR_COPY(X, Y)
#define FLENS_BLASLOG_ERROR_MCOPY(TRANS, X, Y)
#define FLENS_BLASLOG_ERROR_AXPY(ALPHA, X, Y)
#define FLENS_BLASLOG_ERROR_MAXPY(TRANS, ALPHA, X, Y)

#define FLENS_BLASLOG_RESIZE_VECTOR(X, n)
#define FLENS_BLASLOG_RESIZE_MATRIX(B, numRows, numCols)
#define FLENS_BLASLOG_RESIZE_GBMATRIX(B, numRows, numCols, numSubDiags, numSuperDiags)
#define FLENS_BLASLOG_RESIZE_HBMATRIX(B, dim, numOffDiags)
#define FLENS_BLASLOG_RESIZE_TBMATRIX(B, dim, numOffDiags)
#define FLENS_BLASLOG_RESIZE_SBMATRIX(B, dim, numOffDiags)

#define FLENS_BLASLOG_BEGIN_COPY(X, Y)
#define FLENS_BLASLOG_BEGIN_RESIDUAL(b, A, x, y);
#define FLENS_BLASLOG_BEGIN_MCOPY(TRANS, X, Y)
#define FLENS_BLASLOG_BEGIN_MCOTR(TRANS, X)
#define FLENS_BLASLOG_BEGIN_AXPY(ALPHA, X, Y)
#define FLENS_BLASLOG_BEGIN_MAXPY(TRANS, ALPHA, X, Y)
#define FLENS_BLASLOG_BEGIN_RAXPY(ALPHA, X, Y)
#define FLENS_BLASLOG_BEGIN_MRAXPY(TRANS, ALPHA, X, Y)
#define FLENS_BLASLOG_BEGIN_SCAL(ALPHA, X)
#define FLENS_BLASLOG_BEGIN_RSCAL(ALPHA, X)
#define FLENS_BLASLOG_BEGIN_DOT(X, Y)
#define FLENS_BLASLOG_BEGIN_DOTU(X, Y)

#define FLENS_BLASLOG_BEGIN_GEMV(TRANS, ALPHA, A, X, BETA, Y)
#define FLENS_BLASLOG_BEGIN_GBMV(TRANS, ALPHA, A, X, BETA, Y)
#define FLENS_BLASLOG_BEGIN_HBMV(ALPHA, A, X, BETA, Y)
#define FLENS_BLASLOG_BEGIN_HEMV(ALPHA, A, X, BETA, Y)
#define FLENS_BLASLOG_BEGIN_HPMV(ALPHA, A, X, BETA, Y)
#define FLENS_BLASLOG_BEGIN_SBMV(ALPHA, A, X, BETA, Y)
#define FLENS_BLASLOG_BEGIN_SPMV(ALPHA, A, X, BETA, Y)
#define FLENS_BLASLOG_BEGIN_SYMV(ALPHA, A, X, BETA, Y)
#define FLENS_BLASLOG_BEGIN_TBMV(TRANS, A, X)
#define FLENS_BLASLOG_BEGIN_TRMV(TRANS, A, X)
#define FLENS_BLASLOG_BEGIN_TPMV(TRANS, A, X)

#define FLENS_BLASLOG_BEGIN_GBMM(transA, transB, alpha, A, B, beta, C)
#define FLENS_BLASLOG_BEGIN_GEMM(transA, transB, alpha, A, B, beta, C)
#define FLENS_BLASLOG_BEGIN_TRMM(side, transA, alpha, A, B)
#define FLENS_BLASLOG_BEGIN_SYMM(side, alpha, A, B, beta, C)
