/*
 *   Copyright (c) 2009, Michael Lehn
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

#ifndef FLENS_DEBUG_AUXILIARY_VERBALIZECLOSURE_TCC
#define FLENS_DEBUG_AUXILIARY_VERBALIZECLOSURE_TCC 1

#include <sstream>
#include <flens/blas/operators/operators.h>
#include <flens/debug/auxiliary/operation.h>
#include <flens/debug/auxiliary/typeid.h>
#include <flens/debug/auxiliary/verbalizeclosure.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace verbose {

template <typename T>
typename RestrictTo<!IsScalar<T>::value, std::string>::Type
verbalizeClosure(VariablePool &variablePool, const T &x)
{
    return variablePool.name(x);
}

template <typename T>
typename RestrictTo<IsScalar<T>::value, std::string>::Type
verbalizeClosure(VariablePool &, const T &x)
{
    std::stringstream sstream;
    sstream << x;
    return sstream.str();
}

template <typename T>
std::string
verbalizeClosure(VariablePool &variablePool, const ScalarValue<T> &x)
{
    std::stringstream sstream;
    sstream << x.value();
    return sstream.str();
}

template <typename I>
std::string
verbalizeClosure(VariablePool &variablePool, const Matrix<I> &x)
{
    return verbalizeClosure(variablePool, x.impl());
}

template <typename I>
std::string
verbalizeClosure(VariablePool &variablePool, const Vector<I> &x)
{
    return verbalizeClosure(variablePool, x.impl());
}

template <typename Op, typename L, typename R>
std::string
verbalizeClosure(VariablePool &variablePool, const VectorClosure<Op, L, R> &x)
{
    return operation<Op>(verbalizeClosure(variablePool, x.left()),
                         verbalizeClosure(variablePool, x.right()));
}

template <typename Op, typename L, typename R>
std::string
verbalizeClosure(VariablePool &variablePool, const MatrixClosure<Op, L, R> &x)
{
    return operation<Op>(verbalizeClosure(variablePool, x.left()),
                         verbalizeClosure(variablePool, x.right()));
}

} } // namespace verbose, namespace flens

#endif // FLENS_DEBUG_AUXILIARY_VERBALIZECLOSURE_TCC
