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

#ifndef FLENS_IO_BANDSTORAGE_SAVE_TCC
#define FLENS_IO_BANDSTORAGE_SAVE_TCC 1


#include <cxxstd/fstream.h>
#include <cxxstd/iomanip.h>

#include <flens/auxiliary/iscomplex.h>
#include <flens/auxiliary/isinteger.h>
#include <flens/io/bandstorage/save.h>
#include <flens/typedefs.h>

namespace flens {

template <typename FS>
bool
save(std::string filename, const GbMatrix<FS> &A)
{

    typedef typename FS::IndexType   IndexType;
    typedef typename FS::ElementType ElementType;

    std::ofstream ofs( filename.c_str(), std::ios::binary );

    if (ofs.is_open()==false) {
        return false;
    }

    IndexType numRows = A.numRows(), numCols = A.numCols();
    IndexType firstRow = A.firstRow(), firstCol = A.firstCol();
    IndexType numSubDiags = A.numSubDiags(), numSuperDiags = A.numSuperDiags();

    ofs.write(reinterpret_cast<char*>(&numRows), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&numCols), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&numSubDiags), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&numSuperDiags), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&firstRow), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&firstCol), sizeof(IndexType));

    for (IndexType i=-A.numSubDiags(); i<=A.numSuperDiags(); ++i) {
        const auto Diag = A.viewDiag(i);
        for (IndexType j=Diag.firstIndex(); j<=Diag.lastIndex(); ++j) {
            ofs.write(reinterpret_cast<const char*>(&(Diag(j))),
                      sizeof(ElementType));
        }
    }
    ofs.close();
    return true;
}

template <typename FS>
bool
save(std::string filename, const HbMatrix<FS> &A)
{

    typedef typename FS::IndexType   IndexType;
    typedef typename FS::ElementType ElementType;

    std::ofstream ofs( filename.c_str(), std::ios::binary );

    if (ofs.is_open()==false) {
        return false;
    }

    IndexType dim = A.dim();
    IndexType firstIndex = A.firstIndex();
    IndexType numOffDiags = A.numOffDiags();

    ofs.write(reinterpret_cast<char*>(&dim), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&numOffDiags), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&firstIndex), sizeof(IndexType));

    for (IndexType i=-A.numOffDiags(); i<=0; ++i) {
        if (A.upLo()==Lower) {
            const auto Diag = A.viewDiag(i);
            for (IndexType j=Diag.firstIndex(); j<=Diag.lastIndex(); ++j) {
                ofs.write(reinterpret_cast<const char*>(&(Diag(j))),
                          sizeof(ElementType));
            }
        } else {
            const auto Diag = A.viewDiag(-i);
            for (IndexType j=Diag.firstIndex(); j<=Diag.lastIndex(); ++j) {
                ElementType alpha = conjugate(Diag(j));
                ofs.write(reinterpret_cast<const char*>(&alpha),
                          sizeof(ElementType));
            }
        }
    }
    ofs.close();
    return true;
}

template <typename FS>
bool
save(std::string filename, const SbMatrix<FS> &A)
{

    typedef typename FS::IndexType   IndexType;
    typedef typename FS::ElementType ElementType;

    std::ofstream ofs( filename.c_str(), std::ios::binary );

    if (ofs.is_open()==false) {
        return false;
    }

    IndexType dim = A.dim();
    IndexType firstIndex = A.firstIndex();
    IndexType numOffDiags = A.numOffDiags();

    ofs.write( reinterpret_cast<char*>(&dim), sizeof(IndexType) );
    ofs.write( reinterpret_cast<char*>(&numOffDiags), sizeof(IndexType) );
    ofs.write( reinterpret_cast<char*>(&firstIndex), sizeof(IndexType) );

    for (IndexType i=-A.numOffDiags(); i<=0; ++i) {
        if (A.upLo()==Lower) {
            const auto Diag = A.viewDiag(i);

            for (IndexType j=Diag.firstIndex(); j<=Diag.lastIndex(); ++j) {
                ofs.write(reinterpret_cast<const char*>(&(Diag(j))),
                          sizeof(ElementType));
            }
        } else {
            const auto Diag = A.viewDiag(-i);

            for (IndexType j=Diag.firstIndex(); j<=Diag.lastIndex(); ++j) {
                 ofs.write(reinterpret_cast<const char*>(&(Diag(j))),
                           sizeof(ElementType));
            }
        }
    }

    ofs.close();
    return true;
}

template <typename FS>
bool
save(std::string filename, const TbMatrix<FS> &A)
{
    typedef typename FS::IndexType   IndexType;
    typedef typename FS::ElementType ElementType;

    std::ofstream ofs( filename.c_str(), std::ios::binary );

    if (ofs.is_open()==false) {
        return false;
    }

    IndexType dim = A.dim();
    IndexType numOffDiags = A.numOffDiags();
    IndexType firstIndex = A.firstIndex();

    StorageUpLo  upLo = A.upLo();
    Diag         diag = A.diag();

    ofs.write(reinterpret_cast<char*>(&dim), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&numOffDiags), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&firstIndex), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&upLo), sizeof(StorageUpLo));
    ofs.write(reinterpret_cast<char*>(&diag), sizeof(Diag));

    if (upLo == Lower) {
        for (IndexType i=-A.numOffDiags(); i <= 0; ++i){
            const auto Diag = A.viewDiag(i);
            for (IndexType j=Diag.firstIndex(); j<=Diag.lastIndex(); ++j) {
                if (diag==NonUnit || i != 0) {
                    ofs.write(reinterpret_cast<const char*>(&(Diag(j))),
                              sizeof(ElementType) );
                }
            }
        }
    } else {
        for (IndexType i=0; i<=A.numOffDiags(); ++i){
            const auto Diag = A.viewDiag(i);
            for (IndexType j=Diag.firstIndex(); j<=Diag.lastIndex(); ++j) {
                if (diag==NonUnit || i != 0) {
                    ofs.write(reinterpret_cast<const char*>(&(Diag(j))),
                              sizeof(ElementType) );
                }
            }
        }
    }

    ofs.close();
    return true;
}

} // namespace flens

#endif // FLENS_IO_BANDSTORAGE_SAVE_TCC

