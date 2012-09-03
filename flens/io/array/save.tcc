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

#ifndef FLENS_IO_ARRAY_SAVE_TCC
#define FLENS_IO_ARRAY_SAVE_TCC 1

#include <fstream>
#include <flens/io/array/save.h>

namespace flens {

template <typename A>
bool
save(std::string filename, const DenseVector<A> &x)
{

    typedef typename A::IndexType   IndexType;
    typedef typename A::ElementType ElementType;

    std::ofstream ofs( filename.c_str(), std::ios::binary );

    if (ofs.is_open() == false)
        return false;

    IndexType length     = x.length();
    IndexType firstIndex = x.firstIndex();

    ofs.write( reinterpret_cast<char*>(&length), sizeof(IndexType) );
    ofs.write( reinterpret_cast<char*>(&firstIndex), sizeof(IndexType) );

    for (IndexType i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        ofs.write( reinterpret_cast<const char*>(&(x(i))), sizeof(ElementType) );
    }

    ofs.close();
    return true;
}

//-- forwarding ---------------------------------------------------------------

template <typename V>
typename RestrictTo<IsVector<V>::value,
                    bool>::Type
save(std::string filename, const V &&x)
{
    return save(filename, x);
}

} // namespace flens

#endif // FLENS_IO_ARRAY_SAVE_TCC
