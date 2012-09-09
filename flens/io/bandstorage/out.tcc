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

#ifndef FLENS_IO_BANDSTORAGE_OUT_TCC
#define FLENS_IO_BANDSTORAGE_OUT_TCC

#include <cxxblas/typedefs.h>
#include <flens/io/bandstorage/out.h>

namespace flens {

template <typename FS>
std::ostream &
operator<<(std::ostream &out, const GbMatrix<FS> &A)
{
    typedef typename GbMatrix<FS>::IndexType IndexType;

#   ifdef FLENS_IO_WITH_RANGES
    IndexType defaultIndexBase = FS::defaultIndexBase;

    if ((A.firstRow()==defaultIndexBase)
     && (A.firstCol()==defaultIndexBase))
    {
        out << "[" << A.numRows() << ", " << A.numCols()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstRow() << ".." << A.lastRow()
            << ", " << A.firstCol() << ".." << A.lastCol()
            << "]" << std::endl;
    }
#   endif // FLENS_IO_WITH_RANGES

    out << std::endl;
    out.setf(std::ios::fixed|std::ios::right);
    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            if(IsNotComplex<typename FS::ElementType>::value)
                out.width(11);
            else
                out.width(22);
            
            if ((i-j<=A.numSubDiags()) && (j-i<=A.numSuperDiags())) {
                out << A(i,j);
            } else {
                out << " ";
            }
        }
        out << std::endl;
    }
    return out;
}

template <typename FS>
std::ostream &
operator<<(std::ostream &out, const HbMatrix<FS> &A)
{
    typedef typename HbMatrix<FS>::ElementType  ElementType;
    typedef typename HbMatrix<FS>::IndexType    IndexType;
#   ifdef FLENS_IO_WITH_RANGES
    IndexType defaultIndexBase = FS::defaultIndexBase;

    if ((A.firstIndex()==defaultIndexBase)
     && (A.firstIndex()==defaultIndexBase))
    {
        out << "[" << A.dim() << ", " << A.dim()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstIndex() << ".." << A.lastIndex()
            << ", " << A.firstIndex() << ".." << A.lastIndex()
            << "]" << std::endl;
    }
#   endif // FLENS_IO_WITH_RANGES

    out << std::endl;
    out.setf(std::ios::fixed|std::ios::right);
    for (IndexType i=A.firstIndex(); i<=A.lastIndex(); ++i) {
        for (IndexType j=A.firstIndex(); j<=A.lastIndex(); ++j) {
             if(IsNotComplex<typename FS::ElementType>::value)
                out.width(11);
            else
                out.width(22);
            
            if ((i-j<=A.numOffDiags()) && (j-i<=A.numOffDiags())) {
                if (i==j) {
                    out << ElementType(cxxblas::real(A(i,j)));
                }
                if (i<j) {
                    out << ((A.upLo()==cxxblas::Upper)
                            ? A(i,j)
                            : cxxblas::conjugate(A(j,i)));
                }
                if (i>j) {
                    out << ((A.upLo()==cxxblas::Upper)
                            ? cxxblas::conjugate(A(j,i))
                            : A(i,j));
                }
            } else {
                out << " ";
            }
        }
        out << std::endl;
    }
    return out;
}

template <typename FS>
std::ostream &
operator<<(std::ostream &out, const SbMatrix<FS> &A)
{
    typedef typename SbMatrix<FS>::IndexType IndexType;

#   ifdef FLENS_IO_WITH_RANGES
    IndexType defaultIndexBase = FS::defaultIndexBase;

    if ((A.firstIndex()==defaultIndexBase)
     && (A.firstIndex()==defaultIndexBase))
    {
        out << "[" << A.dim() << ", " << A.dim()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstIndex() << ".." << A.lastIndex()
            << ", " << A.firstIndex() << ".." << A.lastIndex()
            << "]" << std::endl;
    }
#   endif // FLENS_IO_WITH_RANGES

    out << std::endl;
    out.setf(std::ios::fixed|std::ios::right);
    for (IndexType i=A.firstIndex(); i<=A.lastIndex(); ++i) {
        for (IndexType j=A.firstIndex(); j<=A.lastIndex(); ++j) {
            if(IsNotComplex<typename FS::ElementType>::value)
                out.width(11);
            else
                out.width(22);
            
            if ((i-j<=A.numOffDiags()) && (j-i<=A.numOffDiags())) {
                out << ((A.upLo()==cxxblas::Upper)
                        ? A(std::min(i,j), std::max(i,j))
                        : A(std::max(i,j), std::min(i,j)));
            } else {
                out << " ";
            }
        }
        out << std::endl;
    }
    return out;
}

template <typename FS>
std::ostream &
operator<<(std::ostream &out, const TbMatrix<FS> &A)
{   
    typedef typename TbMatrix<FS>::IndexType    IndexType;
    typedef typename TbMatrix<FS>::ElementType  ElementType;

#   ifdef FLENS_IO_WITH_RANGES
    IndexType defaultIndexBase = FS::defaultIndexBase;

    if ((A.firstIndex()==defaultIndexBase)
     && (A.firstIndex()==defaultIndexBase))
    {
        out << "[" << A.dim() << ", " << A.dim()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstIndex() << ".." << A.lastIndex()
            << ", " << A.firstIndex() << ".." << A.lastIndex()
            << "]" << std::endl;
    }
#   endif // FLENS_IO_WITH_RANGES

    out << std::endl;
    out.setf(std::ios::fixed|std::ios::right);
    for (IndexType i=A.firstIndex(); i<=A.lastIndex(); ++i) {
        for (IndexType j=A.firstIndex(); j<=A.lastIndex(); ++j) {
            if(IsNotComplex<typename FS::ElementType>::value)
                out.width(11);
            else
                out.width(22);
            
            if (i==j) {
                (A.diag()==cxxblas::Unit) ? out << ElementType(1)
                                          : out << A(i,j);
            } else {
                if (((i>j) && (A.upLo()==cxxblas::Lower) && (i-j<=A.numOffDiags()))
                 || ((i<j) && (A.upLo()==cxxblas::Upper) && (j-i<=A.numOffDiags()))) {
                    out << A(i,j);
                } else {
                    out << " ";
                }
            }
            out << " ";
        }
        out << std::endl;
    }
    return out;
}

} // namespace flens

#endif // FLENS_IO_BANDSTORAGE_OUT_TCC
