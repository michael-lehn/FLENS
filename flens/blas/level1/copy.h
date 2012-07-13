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

#ifndef FLENS_BLAS_LEVEL1_COPY_H
#define FLENS_BLAS_LEVEL1_COPY_H 1

#include <cxxblas/cxxblas.h>
#include <flens/auxiliary/auxiliary.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//-- copy
template <typename VX, typename VY>
    void
    copy(const DenseVector<VX> &x, DenseVector<VY> &y);

//-- gbcopy
template <typename MA, typename MB>
    void
    copy(Transpose trans, const GbMatrix<MA> &A, GbMatrix<MB> &B);

//-- gecopy
template <typename MA, typename MB>
    void
    copy(Transpose trans, const GeMatrix<MA> &A, GeMatrix<MB> &B);

//-- hbcopy
template <typename MA, typename MB>
    void
    copy(Transpose trans, const HbMatrix<MA> &A, HbMatrix<MB> &B);

//-- hpcopy
template <typename MA, typename MB>
    void
    copy(Transpose trans, const HpMatrix<MA> &A, HpMatrix<MB> &B);
    
//-- tbcopy
template <typename MA, typename MB>
    void
    copy(Transpose trans, const TbMatrix<MA> &A, TbMatrix<MB> &B);    

//-- trcopy
template <typename MA, typename MB>
    void
    copy(Transpose trans, const TrMatrix<MA> &A, TrMatrix<MB> &B);

//-- tpcopy
template <typename MA, typename MB>
    void
    copy(Transpose trans, const TpMatrix<MA> &A, TpMatrix<MB> &B);

//-- sbcopy
template <typename MA, typename MB>
    void
    copy(const SbMatrix<MA> &A, SbMatrix<MB> &B);

//-- spcopy
template <typename MA, typename MB>
    void
    copy(const SpMatrix<MA> &A, SpMatrix<MB> &B);

//-- sycopy
template <typename MA, typename MB>
    void
    copy(const SyMatrix<MA> &A, SyMatrix<MB> &B);

//-- extensions ----------------------------------------------------------------

//-- copy: HbMatrix -> GbMatrix
template <typename MA, typename MB>
    void
    copy(const HbMatrix<MA> &A, GbMatrix<MB> &B);  

//-- copy: TbMatrix -> GbMatrix
template <typename MA, typename MB>
    void
    copy(Transpose trans, const TbMatrix<MA> &A, GbMatrix<MB> &B);  
    
//-- copy: TrMatrix -> GeMatrix
template <typename MA, typename MB>
    void
    copy(Transpose trans, const TrMatrix<MA> &A, GeMatrix<MB> &B);

//-- copy: GeMatrix -> TrMatrix
template <typename MA, typename MB>
    void
    copy(Transpose trans, const GeMatrix<MA> &A, TrMatrix<MB> &B);
   
//-- copy: SbMatrix -> GbMatrix
template <typename MA, typename MB>
    void
    copy(const SbMatrix<MA> &A, GbMatrix<MB> &B);

//-- copy: SyMatrix -> GeMatrix
template <typename MA, typename MB>
    void
    copy(const SyMatrix<MA> &A, GeMatrix<MB> &B);

//-- copy: GeMatrix -> SyMatrix
template <typename MA, typename MB>
    void
    copy(const GeMatrix<MA> &A, SyMatrix<MB> &B);

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_COPY_H
