/*
 *   Copyright (c) 2007, Michael Lehn
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

#ifndef FLENS_BLAS_CLOSURES_ASSIGN_TCC
#define FLENS_BLAS_CLOSURES_ASSIGN_TCC 1

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

#include <flens/blas/blas.h>
#include <flens/matrixtypes/matrix.h>
#include <flens/vectortypes/vector.h>

namespace flens {

//-- vector closures
template <typename VX, typename VY>
void
assign(const Vector<VX> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_ASSIGNMENT(x, y);

    blas::copy(x.impl(), y.impl());

    FLENS_BLASLOG_END;
}

template <typename VX, typename VY>
void
plusAssign(const Vector<VX> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_PLUSASSIGNMENT(x, y);

    typedef typename Vector<VX>::Impl::ElementType  T;
    const T  One(1);

    blas::axpy(One, x.impl(), y.impl());

    FLENS_BLASLOG_END;
}

template <typename VX, typename VY>
void
minusAssign(const Vector<VX> &x, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_MINUSASSIGNMENT(x, y);

    typedef typename Vector<VX>::Impl::ElementType  T;
    const T  MinusOne(-1);

    blas::axpy(MinusOne, x.impl(), y.impl());

    FLENS_BLASLOG_END;
}

//-- matrix closures
template <typename MA, typename MB>
void
assign(const Matrix<MA> &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_BEGIN_ASSIGNMENT(A, B);

    blas::copy(NoTrans, A.impl(), B.impl());

    FLENS_BLASLOG_END;
}

template <typename MA, typename MB>
void
plusAssign(const Matrix<MA> &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_BEGIN_PLUSASSIGNMENT(A, B);

    typedef typename Matrix<MA>::Impl::ElementType  T;
    const T  One(1);

    blas::axpy(NoTrans, One, A.impl(), B.impl());

    FLENS_BLASLOG_END;
}

template <typename MA, typename MB>
void
minusAssign(const Matrix<MA> &A, Matrix<MB> &B)
{
    FLENS_BLASLOG_BEGIN_MINUSASSIGNMENT(A, B);

    typedef typename Matrix<MA>::Impl::ElementType  T;
    const T  MinusOne(-1);

    blas::axpy(NoTrans, MinusOne, A.impl(), B.impl());

    FLENS_BLASLOG_END;
}

} // namespace flens

#endif // FLENS_BLAS_CLOSURES_ASSIGN_TCC
