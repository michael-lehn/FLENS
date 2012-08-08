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

#ifndef FLENS_DEBUG_AUXILIARY_OPERATION_TCC
#define FLENS_DEBUG_AUXILIARY_OPERATION_TCC 1

#include <flens/auxiliary/issame.h>
#include <flens/blas/operators/operators.h>

namespace flens { namespace verbose {

template <typename Op>
std::string
operation(const std::string &leftOperand, const std::string &rightOperand)
{
    std::ostringstream s;

    if (IsSame<Op, OpAdd>::value)
    {
        s << "("<< leftOperand << " + " << rightOperand << ")";
        return s.str();
    }
    if (IsSame<Op, OpSub>::value)
    {
        s << "("<< leftOperand << " - " << rightOperand << ")";
        return s.str();
    }
    if (IsSame<Op, OpMult>::value)
    {
        s << "(" << leftOperand << " * " << rightOperand << ")";
        return s.str();
    }
    if (IsSame<Op, OpDiv>::value)
    {
        s << "(" << leftOperand << "/" << rightOperand <<  ")";
        return s.str();
    }
    if (IsSame<Op, OpTrans>::value)
    {
        s << "(" << leftOperand << ")^T";
        return s.str();
    }
    if (IsSame<Op, OpConj>::value)
    {
        s << "conj(" << leftOperand << ")";
        return s.str();
    }
    s << "unknown_op(" << leftOperand << ", " << rightOperand << ")";
    return s.str();
}

} } // namespace verbose, namespace flens

#endif // FLENS_DEBUG_AUXILIARY_OPERATION_TCC
