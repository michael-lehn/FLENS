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

#define FLENS_BLASLOG_VARDECL(x)    x

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_SETTAG(tag)                                           \
    verbose::ClosureLog::setTag(tag)

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_UNSETTAG                                              \
    verbose::ClosureLog::unsetTag()

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_END                                                   \
    verbose::ClosureLog::closeEntry();

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_ASSIGNMENT(X, Y)                                \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::separator();                                   \
        verbose::ClosureLog::append() << Y << " = " << X << ";";            \
        verbose::ClosureLog::append() << "flens::assign(" << X << ", " << Y \
                                      << ");";                              \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_PLUSASSIGNMENT(X, Y)                            \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::separator();                                   \
        verbose::ClosureLog::append() << Y << " += " << X << ";";           \
        verbose::ClosureLog::append() << "flens::plusAssign("               \
                                      << X << ", " << Y << ");";            \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_MINUSASSIGNMENT(X, Y)                           \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::separator();                                   \
        verbose::ClosureLog::append() << Y << " -= " << X << ";";           \
        verbose::ClosureLog::append() << "flens::minusAssign("              \
                                      << X << ", " << Y << ");";            \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_IDENTICAL(A, B)                                       \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "Note that "                       \
                                      << A << " and " << B                  \
                                      << " are identical.";                 \
        verbose::ClosureLog::closeEntry();                                  \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_TMP_ADD(tmp)                                          \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::variablePool.addTemporary(tmp);                \
        verbose::ClosureLog::closeEntry();                                  \
    }


//------------------------------------------------------------------------------

#define FLENS_BLASLOG_TMP_TRON                                              \
    verbose::ClosureLog::variablePool.tmpTron = true;                       \

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_TMP_TROFF                                             \
    verbose::ClosureLog::variablePool.tmpTron = false;                      \

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_TMP_REMOVE(tmp, X)                                    \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "WARNING: Temporary " << tmp       \
                                      << " was required for " << X;         \
        verbose::ClosureLog::variablePool.removeTemporary(tmp);             \
        verbose::ClosureLog::closeEntry();                                  \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_ERROR_COPY(X, Y)                                      \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "ERROR: Can not evaluate";         \
        verbose::ClosureLog::append() << "flens::blas::copy("               \
                                      << X << ", " << Y  << ");";           \
        verbose::ClosureLog::append() << "Unknown operation in: " << X;     \
        verbose::ClosureLog::stop();                                        \
        ASSERT(0);                                                          \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_ERROR_MCOPY(TRANS, X, Y)                              \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "ERROR: Can not evaluate";         \
        verbose::ClosureLog::append() << "flens::blas::copy(" << TRANS      \
                                      << ", " << X << ", " << Y  << ");";   \
        verbose::ClosureLog::append() << "Unknown operation in: " << X;     \
        verbose::ClosureLog::stop();                                        \
        ASSERT(0);                                                          \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_ERROR_AXPY(ALPHA, X, Y)                               \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "ERROR: Can not evaluate";         \
        verbose::ClosureLog::append() << "flens::blas::axpy("               \
                                      << ALPHA << ", " << X << ", " << Y    \
                                      << ");";                              \
        verbose::ClosureLog::append() << "Unknown operation in: " << X;     \
        verbose::ClosureLog::stop();                                        \
        ASSERT(0);                                                          \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_ERROR_MAXPY(TRANS, ALPHA, X, Y)                       \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "ERROR: Can not evaluate";         \
        verbose::ClosureLog::append() << "flens::blas::axpy("               \
                                      << TRANS << ","                       \
                                      << ALPHA << ", " << X << ", " << Y    \
                                      << ");";                              \
        verbose::ClosureLog::append() << "Unknown operation in: " << X;     \
        verbose::ClosureLog::stop();                                        \
        ASSERT(0);                                                          \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_RESIZE_VECTOR(X, n)                                   \
    if (verbose::ClosureLog::openEntry()) {                                 \
        verbose::ClosureLog::append() << "WARNING: Resizing " << X          \
                                      << " old dim = " << X.length()        \
                                      << " new dim = " << n;                \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_RESIZE_MATRIX(A, m, n)                                \
    if (verbose::ClosureLog::openEntry()) {                                 \
        verbose::ClosureLog::append() << "WARNING: Resizing " << A          \
                 << " old dim = " << A.numRows() << " x " << A.numCols()    \
                 << ", new dim = " << m << " x " << n;                      \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_RESIZE_GBMATRIX(A, m, n, kl, ku)                      \
    if (verbose::ClosureLog::openEntry()) {                                 \
        verbose::ClosureLog::append() << "WARNING: Resizing " << A          \
                 << " old dim = " << A.numRows() << " x " << A.numCols()    \
                 << A.numSubDiags() << " SubDiags, " << A.numSuperDiags()   \
                 << " SuperDiags, "                                         \
                 << ", new dim = " << m << " x " << n << ", "               \
                 << kl << " SubDiags, " << ku << " SuperDiags";             \
    }
//------------------------------------------------------------------------------

#define FLENS_BLASLOG_RESIZE_HBMATRIX(A, n, k)                              \
    if (verbose::ClosureLog::openEntry()) {                                 \
        verbose::ClosureLog::append() << "WARNING: Resizing " << A          \
                 << " old dim = " << A.dim() << " x " << A.dim()            \
                 << A.numOffDiags() << " OffDiags, "                        \
                 << ", new dim = " << n << " x " << n << ", "               \
                 << k << " OffDiags ";                                      \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_RESIZE_TBMATRIX(A, n, k)                              \
    if (verbose::ClosureLog::openEntry()) {                                 \
        verbose::ClosureLog::append() << "WARNING: Resizing " << A          \
                 << " old dim = " << A.dim() << " x " << A.dim()            \
                 << A.numOffDiags() << " OffDiags, "                        \
                 << ", new dim = " << n << " x " << n << ", "               \
                 << k << " OffDiags ";                                      \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_RESIZE_SBMATRIX(A, n, k)                              \
    if (verbose::ClosureLog::openEntry()) {                                 \
        verbose::ClosureLog::append() << "WARNING: Resizing " << A          \
                 << " old dim = " << A.dim() << " x " << A.dim()            \
                 << A.numOffDiags() << " OffDiags, "                        \
                 << ", new dim = " << n << " x " << n << ", "               \
                 << k << " OffDiags ";                                      \
    }

//------------------------------------------------------------------------------
#define FLENS_BLASLOG_BEGIN_COPY(CONJ, X, Y)                                \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::copy("               \
                                      << CONJ << ", "                       \
                                      << X << ", " << Y << ");";            \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_RESIDUAL(b, A, x, y)                            \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "residual("                        \
            << b << ", " << A << ", " << x << ", " << y << ");";            \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_MCOPY(TRANS, X, Y)                              \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::copy(" << TRANS      \
                                      << ", " << X << ", " << Y << ");";    \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_MCOTR(TRANS, X)                                 \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::cotr(" << TRANS      \
                                      << ", " << X << ");";                 \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_AXPY(ALPHA, X, Y)                               \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::axpy("               \
                                      << ALPHA << ", " << X << ", " << Y    \
                                      << ");";                              \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_AXPBY(ALPHA, X, BETA, Y)                        \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::axpby("              \
                                      << ALPHA << ", " << X << ", " << BETA \
                                      << ", " << Y << ");";                 \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_MAXPY(TRANS, ALPHA, X, Y)                       \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::axpy(" << TRANS      \
                                      << ", " << ALPHA << ", " << X << ", " \
                                      << Y  << ");";                        \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_MAXPBY(TRANS, ALPHA, X, BETA, Y)                \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::axpby(" << TRANS     \
                                      << ", " << ALPHA << ", " << X << ", " \
                                      << BETA << ", " << Y  << ");";        \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_RAXPY(ALPHA, X, Y)                              \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::raxpy("              \
                                      << ALPHA << ", " << X << ", " << Y    \
                                      << ");";                              \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_MRAXPY(TRANS, ALPHA, X, Y)                      \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::raxpy("              \
                                      << TRANS << ", "                      \
                                      << ALPHA << ", " << X << ", " << Y    \
                                      << ");";                              \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_SCAL(ALPHA, X)                                  \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::scal("               \
                                      << ALPHA << ", " << X << ");";        \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_RSCAL(ALPHA, X)                                 \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::rscal("              \
                                      << ALPHA << ", " << X << ");";        \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_DOT(X, Y)                                       \
   if (verbose::ClosureLog::createEntry()) {                                \
       verbose::ClosureLog::separator();                                    \
       verbose::ClosureLog::append() << "flens::blas::dot("                 \
                                     << X << ", " << Y << ");";             \
}

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_DOTU(X, Y)                                      \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::separator();                                   \
        verbose::ClosureLog::append() << "flens::blas::dotu("               \
                                      << X << ", " << Y << ");";            \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_DOTU(X, Y)                                      \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::separator();                                   \
        verbose::ClosureLog::append() << "flens::blas::dotu("               \
                                      << X << ", " << Y << ");";            \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_GEMV(TRANS, ALPHA, A, X, BETA, Y)               \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mv("                 \
                                      << TRANS << ", " << ALPHA << ", "     \
                                      << A << ", " << X << ", "             \
                                      << BETA << ", "                       \
                                      << Y << ");";                         \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_GBMV(TRANS, ALPHA, A, X, BETA, Y)               \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mv("                 \
                                      << TRANS << ", " << ALPHA << ", "     \
                                      << A << ", " << X << ", "             \
                                      << BETA << ", "                       \
                                      << Y << ");";                         \
    }
//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_HBMV(ALPHA, A, X, BETA, Y)                      \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mv("                 \
                                      << ", " << ALPHA << ", "              \
                                      << A << ", " << X << ", "             \
                                      << BETA << ", "                       \
                                      << Y << ");";                         \
    }
//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_HEMV(ALPHA, A, X, BETA, Y)                      \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mv("                 \
                                      << ", " << ALPHA << ", "              \
                                      << A << ", " << X << ", "             \
                                      << BETA << ", "                       \
                                      << Y << ");";                         \
    }
//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_HPMV(ALPHA, A, X, BETA, Y)                      \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mv("                 \
                                      << ", " << ALPHA << ", "              \
                                      << A << ", " << X << ", "             \
                                      << BETA << ", "                       \
                                      << Y << ");";                         \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_SBMV(ALPHA, A, X, BETA, Y)                      \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mv("                 \
                                      << ", " << ALPHA << ", "              \
                                      << A << ", " << X << ", "             \
                                      << BETA << ", "                       \
                                      << Y << ");";                         \
    }
//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_SPMV(ALPHA, A, X, BETA, Y)                      \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mv("                 \
                                      << TRANS << ", "                      \
                                      << A << ", " << X                     \
                                      << ");";                              \
    }
//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_SYMV(ALPHA, A, X, BETA, Y)                      \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mv("                 \
                                      << ", " << ALPHA << ", "              \
                                      << A << ", " << X << ", "             \
                                      << BETA << ", "                       \
                                      << Y << ");";                         \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_TBMV(TRANS, A, X)                               \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mv("                 \
                                      << TRANS << ", "                      \
                                      << A << ", " << X                     \
                                      << ");";                              \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_TPMV(TRANS, A, X)                               \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mv("                 \
                                      << TRANS << ", "                      \
                                      << A << ", " << X                     \
                                      << ");";                              \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_TRMV(TRANS, A, X)                               \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mv("                 \
                                      << TRANS << ", "                      \
                                      << A << ", " << X                     \
                                      << ");";                              \
    }


//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_GEMM(transA, transB, alpha, A, B, beta, C)      \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mm("                 \
                                      << transA << ", " << transB << ", "   \
                                      << alpha << ", "                      \
                                      << A << ", " << B << ", "             \
                                      << beta << ", "                       \
                                      << C << ");";                         \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_GBMM(transA, transB, alpha, A, B, beta, C)      \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mm("                 \
                                      << transA << ", " << transB << ", "   \
                                      << alpha << ", "                      \
                                      << A << ", " << B << ", "             \
                                      << beta << ", "                       \
                                      << C << ");";                         \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_HBMM(side, alpha, A, B, beta, C)                \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mm("                 \
                                      << side << ", "                       \
                                      << alpha << ", "                      \
                                      << A << ", " << B << ", "             \
                                      << beta << ", "                       \
                                      << C << ");";                         \
    }



//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_HEMM(side, alpha, A, B, beta, C)                \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mm("                 \
                                      << side << ", "                       \
                                      << alpha << ", "                      \
                                      << A << ", " << B << ", "             \
                                      << beta << ", "                       \
                                      << C << ");";                         \
    }



//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_TBMM(side, transA, alpha, A, B)                 \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mm("                 \
                                      << side << ", " << transA << ", "     \
                                      << alpha << ", "                      \
                                      << A << ", " << B << ");";            \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_TRMM(side, transA, alpha, A, B)                 \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mm("                 \
                                      << side << ", " << transA << ", "     \
                                      << alpha << ", "                      \
                                      << A << ", " << B << ");";            \
    }

//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_SBMM(side, alpha, A, B, beta, C)                \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mm("                 \
                                      << side << ", "                       \
                                      << alpha << ", "                      \
                                      << A << ", " << B << ", "             \
                                      << beta << ", "                       \
                                    << C << ");";                         \
    }



//------------------------------------------------------------------------------

#define FLENS_BLASLOG_BEGIN_SYMM(side, alpha, A, B, beta, C)                \
    if (verbose::ClosureLog::createEntry()) {                               \
        verbose::ClosureLog::append() << "flens::blas::mm("                 \
                                      << side << ", "                       \
                                      << alpha << ", "                      \
                                      << A << ", " << B << ", "             \
                                      << beta << ", "                       \
                                      << C << ");";                         \
    }



