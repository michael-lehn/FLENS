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

#include <cxxblas/drivers/drivers.h>

namespace cxxblas {

const char *
blasImpl()
{
    return BLAS_IMPL;
}

//------------------------------------------------------------------------------

namespace CBLAS {

CBLAS_ORDER
value(StorageOrder order)
{
    if (order==RowMajor) {
        return CblasRowMajor;
    }
    return CblasColMajor;
}

CBLAS_TRANSPOSE
value(Transpose trans)
{
    if (trans==NoTrans) {
        return CblasNoTrans;
    }
    if (trans==Conj) {
        return CblasConjNoTrans;
    }
    if (trans==Trans) {
        return CblasTrans;
    }
    return CblasConjTrans;
}

CBLAS_UPLO
value(StorageUpLo upLo)
{
    if (upLo==Upper) {
        return CblasUpper;
    }
    return CblasLower;
}

CBLAS_SIDE
value(Side side)
{
    if (side==Left) {
        return CblasLeft;
    }
    return CblasRight;
}

CBLAS_DIAG
value(Diag diag)
{
    if (diag==Unit) {
        return CblasUnit;
    }
    return CblasNonUnit;
}

} // namespace CBLAS

} // namespace cxxblas
