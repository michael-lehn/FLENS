/*
 *   Copyright (c) 2013, Klaus Pototzky
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

#ifndef PLAYGROUND_FLENS_SPARSE_SUPERLU_CREATEDENSEMATRIX_TCC
#define PLAYGROUND_FLENS_SPARSE_SUPERLU_CREATEDENSEMATRIX_TCC 1

#ifdef WITH_SUPERLU

namespace flens { namespace superlu {

void
Create_Dense_Matrix(SuperMatrix *A, int m, int n, float *dataA, int ldA,
                    Stype_t stype, Mtype_t mtype)
{
    superlu_float::sCreate_Dense_Matrix(A, m, n, dataA, ldA,
                                        stype, SLU_S, mtype);

};

void
Create_Dense_Matrix(SuperMatrix *A, int m, int n, double *dataA, int ldA,
                    Stype_t stype, Mtype_t mtype)
{
    superlu_double::dCreate_Dense_Matrix(A, m, n, dataA, ldA,
                                         stype, SLU_D, mtype);

};

void
Create_Dense_Matrix(SuperMatrix *A, int m, int n,
                    std::complex<float> *dataA, int ldA,
                    Stype_t stype, Mtype_t mtype)
{
    typedef typename superlu_complex_float::complex  slu_complex_float;

    slu_complex_float *dataA_ = reinterpret_cast<slu_complex_float *>(dataA);

    superlu_complex_float::cCreate_Dense_Matrix(A, m, n, dataA_, ldA,
                                                stype, SLU_C, mtype);
};

void
Create_Dense_Matrix(SuperMatrix *A, int m, int n,
                    std::complex<double> *dataA, int ldA,
                    Stype_t stype, Mtype_t mtype)
{
    typedef typename superlu_complex_double::doublecomplex  slu_complex_double;

    slu_complex_double *dataA_ = reinterpret_cast<slu_complex_double *>(dataA);

    superlu_complex_double::zCreate_Dense_Matrix(A, m, n, dataA_, ldA,
                                                 stype, SLU_D, mtype);

};

} } // namespace superlu, flens

#endif // WITH_SUPERLU

#endif // PLAYGROUND_FLENS_SPARSE_SUPERLU_CREATEDENSEMATRIX_TCC
