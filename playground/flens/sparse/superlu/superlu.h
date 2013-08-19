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

#ifndef PLAYGROUND_FLENS_SPARSE_SUPERLU_H
#define PLAYGROUND_FLENS_SPARSE_SUPERLU_H 1

namespace flens { namespace superlu {
    
typedef int int_t;
#include <slu_Cnames.h>
#include <supermatrix.h>
#include <slu_util.h>
    
// Avoid redefinitions in superlu header files
// by using namespaces
    
namespace superlu_float {
#include <slu_sdefs.h>
} // namespace superlu_float
    
namespace superlu_double {
#include <slu_ddefs.h>
} // namespace superlu_double
    
namespace superlu_complex_float{
#include <slu_cdefs.h>
} // namespace superlu_complex_float
    
namespace superlu_complex_double {
#include <slu_zdefs.h>
} // namespace superlu_complex_double
    
} } // namespace superlu, flens


#include<playground/flens/sparse/superlu/create_compcol_matrix.h>
#include<playground/flens/sparse/superlu/create_dense_matrix.h>
#include<playground/flens/sparse/superlu/gssv.h>
#include<playground/flens/sparse/superlu/sv.h>

#endif // PLAYGROUND_FLENS_SPARSE_SUPERLU_H