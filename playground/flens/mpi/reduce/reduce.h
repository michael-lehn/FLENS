/*
 *   Copyright (c) 2012, Michael Lehn, Klaus Pototzky
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

#ifndef PLAYGROUND_FLENS_MPI_REDUCE_REDUCE_H
#define PLAYGROUND_FLENS_MPI_REDUCE_REDUCE_H 1

#include<flens/auxiliary/auxiliary.h>
#include<flens/matrixtypes/matrixtypes.h>
#include<flens/vectortypes/vectortypes.h>
#include<playground/flens/mpi/types.h>


namespace flens { namespace mpi {

#ifdef WITH_MPI

//--- Max ---------------------------------------------------------------------
template <typename T>
    typename RestrictTo<IsReal<T>::value ||
                        IsInteger<T>::value,
                        T>::Type
    MPI_reduce_max(const T &x, const int root = 0,
                   const MPI::Comm &communicator = MPI::COMM_WORLD);

//--- Min ---------------------------------------------------------------------
template <typename T>
    typename RestrictTo<IsReal<T>::value ||
                        IsInteger<T>::value,
                        T>::Type
    MPI_reduce_min(const T &x, const int root = 0,
                   const MPI::Comm &communicator = MPI::COMM_WORLD);

//--- Sum ---------------------------------------------------------------------
template <typename T>
    typename RestrictTo<MPI_Type<T>::Compatible,
                        T>::Type
    MPI_reduce_sum(const T &x, const int root = 0,
                   const MPI::Comm &communicator = MPI::COMM_WORLD);


template <typename VX, typename VSUM>
    typename RestrictTo<IsDenseVector<VX>::value &&
                        IsDenseVector<VSUM>::value,
                        void>::Type
    MPI_reduce_sum(VX &&x, VSUM &&sum, const int root = 0,
                   const MPI::Comm &communicator = MPI::COMM_WORLD);


template <typename MA, typename MSUM>
    typename RestrictTo<IsGeMatrix<MA>::value &&
                        IsGeMatrix<MSUM>::value,
                        void>::Type
    MPI_reduce_sum(MA &&A, MSUM &&Sum, const int root = 0,
                   const MPI::Comm &communicator = MPI::COMM_WORLD);

#else

//--- Max ---------------------------------------------------------------------
template <typename T>
    typename RestrictTo<IsReal<T>::value ||
                        IsInteger<T>::value,
                        T>::Type
    MPI_reduce_max(const T &x, const int root = 0);

//--- Min ---------------------------------------------------------------------
template <typename T>
    typename RestrictTo<IsReal<T>::value ||
                        IsInteger<T>::value,
                        T>::Type
    MPI_reduce_min(const T &x, const int root = 0);

//--- Sum ---------------------------------------------------------------------
template <typename T>
    T
    MPI_reduce_sum(const T &x, const int root = 0);


template <typename X, typename SUM>
    typename RestrictTo<(IsDenseVector<X>::value &&
                         IsDenseVector<SUM>::value) ||
                        (IsGeMatrix<X>::value &&
                         IsGeMatrix<SUM>::value),
                        void>::Type
    MPI_reduce_sum(X &&x, SUM &&sum, const int root = 0);

#endif // WITH_MPI

} }

#endif // PLAYGROUND_FLENS_MPI_REDUCE_REDUCE_H
