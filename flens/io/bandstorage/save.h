/*
 *   Copyright (c) 2012, Klaus Pototzky
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

#ifndef FLENS_IO_BANDSTORAGE_SAVE_H
#define FLENS_IO_BANDSTORAGE_SAVE_H 1

#include <iostream>
#include <limits>
#include <string>

#include <flens/matrixtypes/general/impl/gbmatrix.h>
#include <flens/matrixtypes/hermitian/impl/hbmatrix.h>
#include <flens/matrixtypes/symmetric/impl/sbmatrix.h>
#include <flens/matrixtypes/triangular/impl/tbmatrix.h>

namespace flens {

template <typename FS>
    bool
    save(std::string filename, const GbMatrix<FS> &A);

template <typename FS>
    bool
    save(std::string filename, const HbMatrix<FS> &A);

template <typename FS>
    bool
    save(std::string filename, const SbMatrix<FS> &A);

template <typename FS>
    bool
    save(std::string filename, const TbMatrix<FS> &A);

//-- forwarding ---------------------------------------------------------------
template <typename MA>
    typename RestrictTo<IsGbMatrix<MA>::value ||
                        IsHbMatrix<MA>::value ||
                        IsSbMatrix<MA>::value ||
                        IsTbMatrix<MA>::value,
                        bool>::Type
    save(std::string filename, const MA &&A);

} // namespace flens

#endif // FLENS_IO_BANDSTORAGE_SAVE_H
