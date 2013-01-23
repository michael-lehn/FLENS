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

#ifndef PLAYGROUND_FLENS_MPI_TYPES_H
#define PLAYGROUND_FLENS_MPI_TYPES_H 1

#ifdef WITH_MPI
#    include "mpi.h"
#endif

namespace flens { namespace mpi {


template <typename T>
struct MPI_Type
{
    static const bool            Compatible = false;
};

#ifdef WITH_MPI
template <>
struct MPI_Type<int>
{
    typedef int                   PrimitiveType; 
    static const bool             Compatible = true;
    static MPI_Datatype           Type() { return MPI_INT; }
    static const int              size       = 1;
};

template <>
struct MPI_Type<long>
{
    typedef long                  PrimitiveType; 
    static const bool             Compatible = true;   
    static MPI_Datatype           Type() { return MPI_LONG; }
    static const int              size       = 1;
};

template <>
struct MPI_Type<unsigned int>
{
    typedef int                   PrimitiveType; 
    static const bool             Compatible = true;
    static MPI_Datatype           Type() { return MPI_UNSIGNED; }
    static const int              size       = 1;
};

template <>
struct MPI_Type<unsigned long>
{
    typedef int                   PrimitiveType; 
    static const bool             Compatible = true;   
    static MPI_Datatype           Type() { return MPI_UNSIGNED_LONG; }
    static const int              size       = 1;
};

template <>
struct MPI_Type<float>
{
    typedef float                 PrimitiveType; 
    static const bool             Compatible = true;   
    static MPI_Datatype           Type() { return MPI_FLOAT; }
    static const int              size       = 1;
}; 

template <>
struct MPI_Type<double>
{
    typedef double                PrimitiveType; 
    static const bool             Compatible = true; 
    static MPI_Datatype           Type() { return MPI_DOUBLE; }
    static const int              size       = 1;
};

template <typename T>
struct MPI_Type<std::complex<T> >
{
    typedef T                     PrimitiveType; 
    static const bool             Compatible = true;
    static MPI_Datatype           Type() { return MPI_Type<T>::Type(); }
    static const int              size       = 2*MPI_Type<T>::size;
};

#endif 


} }

#endif // PLAYGROUND_FLENS_MPI_TYPES_H
