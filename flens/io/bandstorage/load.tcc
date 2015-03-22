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

#ifndef FLENS_IO_BANDSTORAGE_LOAD_TCC
#define FLENS_IO_BANDSTORAGE_LOAD_TCC 1

#include <cxxstd/locale.h>
#include <cxxstd/fstream.h>

#include <flens/io/bandstorage/load.h>
#include <flens/typedefs.h>
#include <ulmblas/cxxblas.h>

namespace flens {

template <typename FS>
bool
load(std::string filename, GbMatrix<FS> &A)
{
    typedef typename FS::IndexType   IndexType;
    typedef typename FS::ElementType ElementType;

    std::ifstream ifs(filename.c_str(), std::ios::binary);

    if (ifs.is_open()==false) {
        return false;
    }

    IndexType numRows, numCols;
    IndexType numSubDiags, numSuperDiags;
    IndexType firstRow, firstCol;

    ifs.read(reinterpret_cast<char*>(&numRows), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&numCols), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&numSubDiags), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&numSuperDiags), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&firstRow), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&firstCol), sizeof(IndexType));

    A.resize(numRows, numCols, numSubDiags, numSuperDiags, firstRow, firstCol);

    for (IndexType i=-A.numSubDiags(); i<=A.numSuperDiags(); ++i) {
        auto Diag = A.viewDiag(i);
        for (IndexType j=Diag.firstIndex(); j<=Diag.lastIndex(); ++j) {
            ifs.read(reinterpret_cast<char*>(&(Diag(j))), sizeof(ElementType));
        }
    }

    ifs.close();
    return true;
}


template <typename FS>
bool
load(std::string filename, HbMatrix<FS> &A)
{
    typedef typename FS::IndexType   IndexType;
    typedef typename FS::ElementType ElementType;

    std::ifstream ifs( filename.c_str(), std::ios::binary );

    if (ifs.is_open()==false) {
        return false;
    }

    IndexType dim = A.dim();
    IndexType numOffDiags = A.numOffDiags();
    IndexType firstIndex = A.firstIndex();

    ifs.read(reinterpret_cast<char*>(&dim), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&numOffDiags), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&firstIndex), sizeof(IndexType));

    A.resize(dim, numOffDiags, firstIndex);

    for (IndexType i=-A.numOffDiags(); i<=0; ++i) {
        if (A.upLo()==Lower) {
            auto Diag = A.viewDiag(i);
            for (IndexType j=Diag.firstIndex(); j<=Diag.lastIndex(); ++j) {
                ifs.read(reinterpret_cast<char*>(&(Diag(j))),
                         sizeof(ElementType));
            }
        } else {
            auto Diag = A.viewDiag(-i);
            for (IndexType j=Diag.firstIndex(); j<=Diag.lastIndex(); ++j) {
                ElementType alpha;
                ifs.read(reinterpret_cast<char*>(&alpha), sizeof(ElementType));
                Diag(j) = conjugate(alpha);
            }
        }
    }

    ifs.close();
    return true;
}

template <typename FS>
bool
load(std::string filename, SbMatrix<FS> &A)
{
    typedef typename FS::IndexType   IndexType;
    typedef typename FS::ElementType ElementType;

    std::ifstream ifs( filename.c_str(), std::ios::binary );

    if (ifs.is_open()==false) {
        return false;
    }

    IndexType dim = A.dim();
    IndexType numOffDiags = A.numOffDiags();
    IndexType firstIndex = A.firstIndex();

    ifs.read(reinterpret_cast<char*>(&dim), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&numOffDiags), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&firstIndex), sizeof(IndexType));

    A.resize(dim, numOffDiags, firstIndex);

    for (IndexType i=-A.numOffDiags(); i<=0; ++i) {
        if (A.upLo()==Lower) {
            auto Diag = A.viewDiag(i);

            for (IndexType j=Diag.firstIndex(); j<=Diag.lastIndex(); ++j) {
                ifs.read(reinterpret_cast<char*>(&(Diag(j))),
                         sizeof(ElementType) );
            }
        } else {
            auto Diag = A.viewDiag(-i);

            for (IndexType j=Diag.firstIndex(); j<=Diag.lastIndex(); ++j) {
                 ifs.read(reinterpret_cast<char*>(&(Diag(j))),
                          sizeof(ElementType) );
            }
        }
    }

    ifs.close();
    return true;
}

template <typename FS>
bool
load(std::string filename, TbMatrix<FS> &A)
{
    typedef typename FS::IndexType   IndexType;
    typedef typename FS::ElementType ElementType;

    std::ifstream ifs( filename.c_str(), std::ios::binary );

    if (ifs.is_open()==false) {
        return false;
    }

    IndexType dim;
    IndexType numOffDiags;
    IndexType firstIndex;
    StorageUpLo  upLo;
    Diag         diag;

    ifs.read(reinterpret_cast<char*>(&dim), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&numOffDiags), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&firstIndex), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&upLo), sizeof(StorageUpLo));
    ifs.read(reinterpret_cast<char*>(&diag), sizeof(Diag));

    ASSERT(upLo==A.upLo());
    ASSERT((diag==A.diag()) || (A.diag()==NonUnit));

    A.resize(dim, numOffDiags, firstIndex);

    if (upLo == Lower) {
        for (IndexType i=-A.numOffDiags(); i <= 0; ++i){
            auto Diag = A.viewDiag(i);
            for (IndexType j=Diag.firstIndex(); j<=Diag.lastIndex(); ++j) {
                if (diag==NonUnit || i != 0) {
                    ifs.read(reinterpret_cast<char*>(&(Diag(j))),
                             sizeof(ElementType) );
                } else if (A.diag()==NonUnit) {
                    Diag(i) = ElementType(i);
                }
            }
        }
    } else {
        for (IndexType i=0; i<=A.numOffDiags(); ++i){
            auto Diag = A.viewDiag(i);
            for (IndexType j=Diag.firstIndex(); j<=Diag.lastIndex(); ++j) {
                if (diag==NonUnit || i != 0) {
                    ifs.read(reinterpret_cast<char*>(&(Diag(j))),
                             sizeof(ElementType) );
                } else if (A.diag()==NonUnit) {
                    Diag(i) = ElementType(i);
                }
            }
        }
    }

    ifs.close();
    return true;
}

//-- forwarding ---------------------------------------------------------------
template <typename MA>
typename RestrictTo<IsGbMatrix<MA>::value ||
                    IsHbMatrix<MA>::value ||
                    IsSbMatrix<MA>::value ||
                    IsTbMatrix<MA>::value,
                    bool>::Type
load(std::string filename, MA &&A)
{
    return load(filename, A);
}

} // namespace flens

#endif // FLENS_IO_BANDSTORAGE_LOAD_TCC

