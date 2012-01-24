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

#ifndef FLENS_BLAS_CLOSURES_RESULT_H
#define FLENS_BLAS_CLOSURES_RESULT_H 1

#include <flens/blas/operators/operators.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens {

//-- General definition --------------------------------------------------------
template <typename A>
struct Result
{
    typedef A  Type;
};

//-- Vectors -------------------------------------------------------------------
template <typename V>
struct Result<Vector<V> >
{
    typedef typename Vector<V>::Impl  Type;
};

//== Define results of vector-vector operations ================================
template <typename Op, typename VL, typename VR>
struct Result<VectorClosure<Op, VL, VR> >
{
    typedef typename Result<VL>::Type  _VL;
    typedef typename Result<VR>::Type  _VR;

    typedef typename Result<VectorClosure<Op, _VL, _VR> >::Type Type;
};

//-- DenseVector-DenseVector-operations
//
// Op = OpAdd
//
template <typename VL, typename VR>
struct Result<VectorClosure<OpAdd, DenseVector<VL>, DenseVector<VR> > >
{
    typedef typename Promotion<VL, VR>::Type V;
    typedef typename DenseVector<V>::NoView  Type;
};
//
// Op = OpSub
//
template <typename VL, typename VR>
struct Result<VectorClosure<OpSub, DenseVector<VL>, DenseVector<VR> > >
{
    typedef typename Promotion<VL, VR>::Type V;
    typedef typename DenseVector<V>::NoView  Type;
};


} // namespace flens

#endif // FLENS_BLAS_CLOSURES_PRUNEVECTORCLOSURE_H
