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

#ifndef PLAYGROUND_FLENS_SPARSE_SUPERLU_CREATECOMPCOLMATRIX_TCC
#define PLAYGROUND_FLENS_SPARSE_SUPERLU_CREATECOMPCOLMATRIX_TCC 1

#ifdef WITH_SUPERLU

namespace flens { namespace superlu {

void
Create_CompCol_Matrix(SuperMatrix *A, int numRows, int numCols,
                      int numNonZeros, float *values,
                      int * rows, int *cols, Stype_t order, Mtype_t shape)
{

    superlu_float::sCreate_CompCol_Matrix(A, numRows, numCols, numNonZeros,
                                          values, rows, cols,
                                          order, SLU_S, shape);
};

void
Create_CompCol_Matrix(SuperMatrix *A, int numRows, int numCols,
                      int numNonZeros, double *values,
                      int * rows, int *cols, Stype_t order, Mtype_t shape)
{

    superlu_double::dCreate_CompCol_Matrix(A, numRows, numCols, numNonZeros,
                                           values, rows, cols,
                                           order, SLU_D, shape);
};

void
Create_CompCol_Matrix(SuperMatrix *A, int numRows, int numCols,
                      int numNonZeros, std::complex<float> *values,
                      int * rows, int *cols, Stype_t order, Mtype_t shape)
{
    typedef typename superlu_complex_float::complex  slu_complex_float;

    slu_complex_float *_value = reinterpret_cast<slu_complex_float *>(values);

    superlu_complex_float::cCreate_CompCol_Matrix(A, numRows, numCols,
                                                  numNonZeros,
                                                  _values, rows, cols,
                                                  order, SLU_C, shape);
};

void
Create_CompCol_Matrix(SuperMatrix *A, int numRows, int numCols,
                      int numNonZeros, std::complex<double> *values,
                      int * rows, int *cols, Stype_t order, Mtype_t shape)
{
    typedef typename superlu_complex_double::doublecomplex  slu_complex_double;

    slu_complex_double *_value = reinterpret_cast<slu_complex_double *>(values);

    superlu_complex_double::zCreate_CompCol_Matrix(A, numRows, numCols,
                                                   numNonZeros,
                                                   _values, rows, cols,
                                                   order, SLU_Z, shape);
};

} } // namespace superlu, flens

#endif // WITH_SUPERLU

#endif // PLAYGROUND_FLENS_SPARSE_SUPERLU_CREATECOMPCOLMATRIX_TCC
