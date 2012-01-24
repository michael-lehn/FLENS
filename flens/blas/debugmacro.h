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


#ifndef FLENS_BLAS_DEBUGMACRO_H
#define FLENS_BLAS_DEBUGMACRO_H 1

#ifndef FLENS_DEBUG_CLOSURES
#   define FLENS_CLOSURELOG_END
#   define FLENS_CLOSURELOG_BEGIN_ASSIGNMENT(Y, X)
#   define FLENS_CLOSURELOG_ADD_TMP(tmp, Y)
#   define FLENS_CLOSURELOG_BEGIN_COPY(X, Y)
#   define FLENS_CLOSURELOG_BEGIN_AXPY(ALPHA, X, Y)
#else
#   define FLENS_CLOSURELOG_END                                     \
    verbose::ClosureLog::closeEntry();

#   define FLENS_CLOSURELOG_BEGIN_ASSIGNMENT(X, Y)                  \
    if (verbose::ClosureLog::createEntry()) {                       \
        verbose::ClosureLog::append() << Y << " = " << X << ";";    \
        verbose::ClosureLog::append() << "flens::assign("             \
                                      << X << ", " << Y             \
                                      << ");";                      \
    }

#   define FLENS_CLOSURELOG_ADD_TMP(tmp, X)                         \
    verbose::ClosureLog::variablePool.addTemporary(tmp);            \
    verbose::ClosureLog::append() << "Warning: Temporary needed "   \
                                  << "for " << X << ".";            \
    verbose::ClosureLog::append();

#   define FLENS_CLOSURELOG_BEGIN_COPY(X, Y)                        \
    if (verbose::ClosureLog::createEntry()) {                       \
        verbose::ClosureLog::append() << "flens::blas::copy("       \
                                      << X << ", " << Y             \
                                      << ");";                      \
    }

#   define FLENS_CLOSURELOG_BEGIN_AXPY(ALPHA, X, Y)                 \
    if (verbose::ClosureLog::createEntry()) {                       \
        verbose::ClosureLog::append() << "flens::blas::axpy("       \
                                      << ALPHA << ", "              \
                                      << X << ", " << Y             \
                                      << ");";                      \
    }
#endif // FLENS_DEBUG_CLOSURES

#endif // FLENS_BLAS_DEBUGMACRO_H
