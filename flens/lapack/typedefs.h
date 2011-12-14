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

#ifndef FLENS_LAPACK_TYPEDEFS_H
#define FLENS_LAPACK_TYPEDEFS_H 1

#include <cassert>
#include <complex>
#include <cxxblas/cxxblas.h>
#include <flens/storage/storage.h>

namespace flens { namespace lapack {


// vector views
template <typename T>
using DenseVectorConstView = DenseVector<ConstArrayView<T> >;

template <typename T>
using DenseVectorView = DenseVector<ArrayView<T> >;


// matrix views
template <typename T>
using GeMatrixConstView = GeMatrix<ConstFullStorageView<T, ColMajor> >;

template <typename T>
using GeMatrixView = GeMatrix<FullStorageView<T, ColMajor> >;


enum Norm {
    OneNorm = 'O',
    InfinityNorm = 'I',
    FrobeniusNorm = 'F',
    MaximumNorm = 'M'
};

enum ProbabilityDistribution {
    Uniform01,
    Uniform_11,
    StandardNormal
};

enum Direction {
    Forward,
    Backward
};

enum StoreVectors {
    ColumnWise,
    RowWise
};

enum MachineParameter {
    Eps,               // relative machine precision
    SafeMin,           // safe minimum, such that reciprocal does not overflow
    Base,              // base of the machine
    Precision,         // Eps*Base
    Mantissa,          // number of (base) digits in the mantissa
    Rounding,          // return 1 when rounding occurs in addition, 0 otherwise
    UnderflowExp,      // minimum exponent before (gradual) underflow
    UnderflowThreshold,// underflow threshold - Base**(UnderflowExp-1)
    OverflowExp,       // largest exponent before overflow
    OverflowThreshold, // overflow threshold - Base**OverflowExp*(1-Eps)
};

namespace BALANCE {

    enum Balance {
        None        = 'N',
        PermuteOnly = 'P',
        ScaleOnly   = 'S',
        Both        = 'B'
    };

}

namespace SENSE {

    enum Sense {
        None                  = 'N',
        EigenvaluesOnly       = 'E',
        InvariantSubspaceOnly = 'V',
        Both                  = 'B'
    };

}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_TYPEDEFS_H
