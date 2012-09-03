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

#ifndef FLENS_IO_FULLSTORAGE_LOAD_TCC
#define FLENS_IO_FULLSTORAGE_LOAD_TCC 1

#include <functional>
#include <locale>
#include <fstream>
#include <sstream>
#include <string>

#include <cxxblas/typedefs.h>
#include <flens/io/fullstorage/load.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens {

template <typename FS>
bool
load(std::string filename, GeMatrix<FS> &A)
{
    typedef typename FS::IndexType    IndexType;
    typedef typename FS::ElementType  ElementType;

    std::ifstream ifs( filename.c_str(), std::ios::binary );

    if (ifs.is_open()==false) {
        return false;
    }

    IndexType numRows, numCols;
    IndexType firstRow, firstCol;

    ifs.read(reinterpret_cast<char*>(&numRows), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&numCols), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&firstRow), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&firstCol), sizeof(IndexType));

    A.resize(numRows, numCols, firstRow, firstCol);

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            ifs.read( reinterpret_cast<char*>(&A(i,j)), sizeof(ElementType) );
        }
    }

    ifs.close();
    return true;
}


template <typename FS>
bool
load(std::string filename, HeMatrix<FS> &A)
{
    typedef typename FS::IndexType    IndexType;
    typedef typename FS::ElementType  ElementType;

    std::ifstream ifs( filename.c_str(), std::ios::binary );

    if (ifs.is_open()==false) {
        return false;
    }

    IndexType dim = A.dim();
    IndexType firstIndex = A.firstRow();

    ifs.read(reinterpret_cast<char*>(&dim), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&firstIndex), sizeof(IndexType));

    A.resize(dim, firstIndex);

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=i; ++j) {
            if (A.upLo()==cxxblas::Lower) {
                ifs.read(reinterpret_cast<char*>(&(A(i,j))),
                         sizeof(ElementType) );
            } else {
                ElementType alpha;
                ifs.read(reinterpret_cast<char*>(&alpha), sizeof(ElementType));
                A(j,i) = cxxblas::conjugate(alpha);
            }
        }
    }

    ifs.close();
    return true;
}


template <typename FS>
bool
load(std::string filename, SyMatrix<FS> &A)
{
    typedef typename FS::IndexType    IndexType;
    typedef typename FS::ElementType  ElementType;

    std::ifstream ifs( filename.c_str(), std::ios::binary );

    if (ifs.is_open()==false) {
        return false;
    }

    IndexType dim = A.dim();
    IndexType firstIndex = A.firstRow();

    ifs.read(reinterpret_cast<char*>(&dim), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&firstIndex), sizeof(IndexType));

    A.resize(dim, firstIndex);

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=i; ++j) {
            if (A.upLo()==cxxblas::Lower) {
                ifs.read(reinterpret_cast<char*>(&(A(i,j))),
                         sizeof(ElementType) );
            } else {
                ifs.read(reinterpret_cast<char*>(&(A(j,i))),
                         sizeof(ElementType) );
            }
        }
    }

    ifs.close();
    return true;
}

template <typename FS>
bool
load(std::string filename, TrMatrix<FS> &A)
{
    typedef typename FS::IndexType    IndexType;
    typedef typename FS::ElementType  ElementType;

    std::ifstream ifs( filename.c_str(), std::ios::binary );

    if (ifs.is_open()==false) {
        return false;
    }

    IndexType   numRows, numCols;
    IndexType   firstRow, firstCol;
    StorageUpLo upLo;
    Diag        diag;

    ifs.read(reinterpret_cast<char*>(&numRows), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&numCols), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&firstRow), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&firstCol), sizeof(IndexType));
    ifs.read(reinterpret_cast<char*>(&upLo), sizeof(StorageUpLo));
    ifs.read(reinterpret_cast<char*>(&diag), sizeof(Diag));

    ASSERT(upLo==A.upLo());
    ASSERT((diag==A.diag()) || (A.diag()==cxxblas::NonUnit));

    A.resize(numRows, numCols, firstRow, firstCol);

    if (upLo==cxxblas::Lower) {
        for (IndexType i=A.firstRow(); i <= A.lastRow(); ++i){
            const IndexType jmax = i-A.firstRow()+A.firstCol();

            for (IndexType j=A.firstCol(); j<jmax; ++j) {
                ifs.read(reinterpret_cast<char*>(&(A(i,j))),
                         sizeof(ElementType) );
            }

            if (diag==cxxblas::NonUnit) {
                ifs.read(reinterpret_cast<char*>(&(A(i,jmax))),
                         sizeof(ElementType) );
            } else if ((A.diag()==cxxblas::NonUnit) && (diag==cxxblas::Unit)) {
                A(i,jmax) = ElementType(1);
            }
        }
    } else {
        for (IndexType i=A.firstRow(); i <= A.lastRow(); ++i){
            const IndexType jmin = i-A.firstRow()+A.firstCol();

            if (diag==cxxblas::NonUnit) {
                ifs.read(reinterpret_cast<char*>(&(A(i,jmin))),
                         sizeof(ElementType) );
            } else if ((A.diag()==cxxblas::NonUnit) && (diag==cxxblas::Unit)) {
                A(i, jmin) = ElementType(1);
            }

            for (IndexType j=jmin+1; j<=A.lastCol(); ++j) {
                ifs.read(reinterpret_cast<char*>(&(A(i,j))),
                         sizeof(ElementType));
            }

        }
    }

    ifs.close();
    return true;
}

template <typename FS>
typename RestrictTo<IsReal<typename FS::ElementType>::value, bool>::Type
loadMatrixMarket(std::string filename, GeMatrix<FS> &A)
{
    using std::string;

    typedef typename FS::IndexType    IndexType;
    typedef typename FS::ElementType  ElementType;

    string line;

    std::ifstream ifs( filename.c_str(), std::ios::in );

    if (ifs.is_open()==false) {
        return false;
    }

    IndexType numRows, numCols;

    std::getline (ifs,line);

    #   ifndef NDEBUG
    // transform line to lower case
    std::transform(line.begin(), line.end(), line.begin(),
             std::bind2nd(std::ptr_fun(&std::tolower<char>), std::locale("")));

    std::stringstream ss (line);
    std::string buf;
    ss >> buf;
    ASSERT(buf=="%%matrixmarket");
    ss >> buf;
    ASSERT(buf=="matrix");
    ss >> buf;
    ASSERT(buf=="array");
    ss >> buf;
    ASSERT((buf=="real") || (buf=="complex"));
    ss >> buf;
    ASSERT(buf=="general");
#   endif

    while (ifs.good() && (line.c_str()[0]=='%')) {
        std::getline (ifs,line);
    }


    std::stringstream sline;
    sline << line;
    sline >> numRows >> numCols;

    A.resize(numRows, numCols);

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            if (ifs.good()) {
                getline(ifs,line);
            } else {
                return false;
            }
            std::stringstream ssline(line);
            ssline >> A(i,j);
        }
    }
    ifs.close();
    return true;
}


template <typename FS>
typename RestrictTo<IsComplex<typename FS::ElementType>::value, bool>::Type
loadMatrixMarket(std::string filename, GeMatrix<FS> &A)
{
    using std::string;

    typedef typename FS::IndexType                            IndexType;
    typedef typename FS::ElementType                          ElementType;
    typedef typename ComplexTrait<ElementType>::PrimitiveType PrimitiveType;

    string line;

    std::ifstream ifs( filename.c_str(), std::ios::in );

    if (ifs.is_open()==false) {
        return false;
    }

    IndexType numRows, numCols;

    std::getline (ifs,line);

    #   ifndef NDEBUG
    // transform line to lower case
    std::transform(line.begin(), line.end(), line.begin(),
                   std::bind2nd(std::ptr_fun(&std::tolower<char>),
                                std::locale("")));

    std::stringstream ss (line);
    std::string buf;
    ss >> buf;
    ASSERT(buf=="%%matrixmarket");
    ss >> buf;
    ASSERT(buf=="matrix");
    ss >> buf;
    ASSERT(buf=="array");
    ss >> buf;
    ASSERT(buf=="complex");
    ss >> buf;
    ASSERT(buf=="general");
#   endif

    while (ifs.good() && (line.c_str()[0]=='%')) {
        std::getline (ifs,line);
    }


    std::stringstream sline;
    sline << line;
    sline >> numRows >> numCols;
    if ((A.numRows()==0) && (A.numCols()==0)) {
        A.resize(numRows, numCols);
    }

#   ifndef NDEBUG
    ASSERT(A.numRows()==numRows);
    ASSERT(A.numCols()==numCols);
#   endif


    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            if (ifs.good()) {
                std::getline (ifs,line);
            } else {
                return false;
            }
            std::stringstream ssline(line);
            PrimitiveType a, b;
            ssline >> a >> b;
            A(i,j) = ElementType(a, b);
        }
    }
    ifs.close();
    return true;
}

template <typename FS>
typename RestrictTo<IsReal<typename FS::ElementType>::value, bool>::Type
loadMatrixMarket(std::string filename, SyMatrix<FS> &A)
{
    using std::string;

    typedef typename FS::IndexType    IndexType;
    typedef typename FS::ElementType  ElementType;

    string line;

    std::ifstream ifs(filename.c_str(), std::ios::in);

    if (ifs.is_open()==false) {
        return false;
    }

    IndexType numRows, numCols;

    std::getline (ifs,line);

#   ifndef NDEBUG
    // transform line to lower case
    std::transform(line.begin(), line.end(), line.begin(),
             std::bind2nd(std::ptr_fun(&std::tolower<char>), std::locale("")));

    std::stringstream ss (line);
    std::string buf;
    ss >> buf;
    ASSERT(buf=="%%matrixmarket");
    ss >> buf;
    ASSERT(buf=="matrix");
    ss >> buf;
    ASSERT(buf=="array");
    ss >> buf;
    ASSERT((buf=="integer") || (buf=="real"));
    ss >> buf;
    ASSERT(buf=="symmetric");
#   endif

    while (ifs.good() && (line.c_str()[0]=='%')) {
        std::getline (ifs,line);
    }

    std::stringstream sline;
    sline << line;
    sline >> numRows >> numCols;

    A.resize(numRows, numCols);

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=i; j<=A.lastCol(); ++j) {
            if (ifs.good()) {
                std::getline (ifs,line);
            } else {
                return false;
            }
            std::stringstream ssline(line);

            if (A.upLo()==Upper) {
                ssline >> A(i,j);
            } else {
                ssline >> A(j,i);
            }
        }
    }
    ifs.close();
    return true;
}


template <typename FS>
typename RestrictTo<IsComplex<typename FS::ElementType>::value, bool>::Type
loadMatrixMarket(std::string filename, SyMatrix<FS> &A)
{
    using std::string;

    typedef typename FS::IndexType                            IndexType;
    typedef typename FS::ElementType                          ElementType;
    typedef typename ComplexTrait<ElementType>::PrimitiveType PrimitiveType;

    string line;

    std::ifstream ifs( filename.c_str(), std::ios::in );

    if (ifs.is_open()==false) {
        return false;
    }

    IndexType numRows, numCols;

    std::getline (ifs,line);

#   ifndef NDEBUG
    // transform line to lower case
    std::transform(line.begin(), line.end(), line.begin(),
                   std::bind2nd(std::ptr_fun(&std::tolower<char>),
                                std::locale("")));

    std::stringstream ss (line);
    std::string buf;
    ss >> buf;
    ASSERT(buf=="%%matrixmarket");
    ss >> buf;
    ASSERT(buf=="matrix");
    ss >> buf;
    ASSERT(buf=="array");
    ss >> buf;
    ASSERT(buf=="complex");
    ss >> buf;
    ASSERT(buf=="symmetric");
#   endif

    while (ifs.good() && (line.c_str()[0]=='%')) {
        std::getline (ifs,line);
    }


    std::stringstream sline;
    sline << line;
    sline >> numRows >> numCols;

    A.resize(numRows, numCols);


    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=i; j<=A.lastCol(); ++j) {
            if (ifs.good()) {
                std::getline (ifs,line);
            } else {
                return false;
            }

            std::stringstream ssline(line);
            PrimitiveType a, b;
            ssline >> a >> b;
            if (A.upLo()==Upper) {
                A(i,j) = ElementType(a, b);
            } else {
                A(j,i) = ElementType(a, b);
            }

        }
    }
    ifs.close();
    return true;
}

template <typename FS>
typename RestrictTo<IsComplex<typename FS::ElementType>::value, bool>::Type
loadMatrixMarket(std::string filename, HeMatrix<FS> &A)
{

    typedef typename FS::IndexType                            IndexType;
    typedef typename FS::ElementType                          ElementType;
    typedef typename ComplexTrait<ElementType>::PrimitiveType PrimitiveType;

    std::string line;
    IndexType numRows, numCols;

    std::ifstream ifs(filename.c_str(), std::ios::in);

    if (ifs.is_open()==false) {
        return false;
    }

    std::getline(ifs,line);

#   ifndef NDEBUG
    // transform line to lower case
    std::transform(line.begin(), line.end(), line.begin(),
             std::bind2nd(std::ptr_fun(&std::tolower<char>), std::locale("")));

    std::stringstream ss (line);
    std::string buf;
    ss >> buf;
    ASSERT(buf=="%%matrixmarket");
    ss >> buf;
    ASSERT(buf=="matrix");
    ss >> buf;
    ASSERT(buf=="array");
    ss >> buf;
    ASSERT(buf=="complex");
    ASSERT( IsComplex<ElementType>::value );
    ss >> buf;
    ASSERT(buf=="hermitian");
#   endif

    while (ifs.good() && (line.c_str()[0]=='%')) {
        std::getline (ifs,line);
    }


    std::stringstream sline;
    sline << line;
    sline >> numRows >> numCols;

    A.resize(numRows, numCols);


    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=i; j<=A.lastCol(); ++j) {
            if (ifs.good()) {
                std::getline(ifs,line);
            } else {
                return false;
            }

            std::stringstream ssline(line);
            PrimitiveType a, b;
            ssline >> a >> b;

            if (A.upLo()==Upper) {
                A(i,j) = ElementType(a, -b);
            } else {
                A(j,i) = ElementType(a, b);
            }

        }
    }
    ifs.close();
    return true;
}


//-- forwarding ---------------------------------------------------------------
template <typename MA>
typename RestrictTo<IsGeMatrix<MA>::value ||
                    IsHeMatrix<MA>::value ||
                    IsSyMatrix<MA>::value ||
                    IsTrMatrix<MA>::value,
                    bool>::Type
load(std::string filename, MA &&A)
{
    return load(filename, A);
}

template <typename MA>
typename RestrictTo<IsGeMatrix<MA>::value ||
                    IsHeMatrix<MA>::value ||
                    IsSyMatrix<MA>::value ||
                    IsTrMatrix<MA>::value,
                    bool>::Type
loadMatrixMarket(std::string filename, MA &&A)
{
    return loadMatrixMarket(filename, A);
}

} // namespace flens

#endif // FLENS_IO_FULLSTORAGE_LOAD_TCC

