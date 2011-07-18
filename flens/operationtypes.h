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

#ifndef FLENS_OPERATIONTYPES_H
#define FLENS_OPERATIONTYPES_H 1

namespace flens {

// == Operation types ==========================================================

//-- used by matvec-closures (and scalar closures) -----------------------------
struct OpAdd {};              // x+y, A+B
struct OpSub {};              // x-y, A-B
struct OpMult {};             // alpha*x, alpha*A, x*y = x^T*y, A*x, A*B

struct OpTrans {};            // transpose(A) = A^T
struct OpConjTrans {};        // conjugateTranspose(A) = A^H

//-- used by scalarclosures-----------------------------------------------------

//-- binary and unary operations
struct OpMinus {};              //  x-y, -x

//-- binary operations
struct OpDiv {};                //  x / y
struct OpMod {};                //  x % y
struct OpMax {};                //  max(x,y)
struct OpMin {};                //  min(x,y)
struct OpPow {};                //  x^y

//-- unary operations
struct OpConj {};               //  conj(x)
struct OpCeil {};               //  ceil(x)
struct OpFloor {};              //  floor(x)
struct OpAbs {};                //  abs(x) = |x|
struct OpSqrt {};               //  sqrt(x)
struct OpExp {};                //  exp(x) = e^x
struct OpExp2 {};               //  exp2(x) = 2^x
struct OpLog {};                //  log(x) = ln(x)
struct OpLog2 {};               //  log2(x) = log_2(x)
struct OpLog10 {};              //  log10(x) = log_10(x)
struct OpCos {};                //  cos(x)
struct OpSin {};                //  sin(x)
struct OpTan {};                //  tan(x)
struct OpCosh {};               //  cosh(x)
struct OpSinh {};               //  sinh(x)
struct OpTanh {};               //  tanh(x)
struct OpArccos {};             //  acos(x) = arccos(x)
struct OpArcsin {};             //  asin(x) = arcsin(x)
struct OpArctan {};             //  atan(x) = arctan(x)
struct OpArcosh {};             //  acosh(x) = arcosh(x)
struct OpArsinh {};             //  asinh(x) = arsinh(x)
struct OpArtanh {};             //  atanh(x) = artanh(x)

} // namespace flens

#endif // FLENS_OPERATIONTYPES_H
