/*
 *   Copyright (c) 2010, Michael Lehn
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

#ifndef FLENS_IO_FULLSTORAGE_OUT_TCC
#define FLENS_IO_FULLSTORAGE_OUT_TCC 1

#include <flens/auxiliary/iscomplex.h>
#include <cxxblas/typedefs.h>
#include <flens/io/fullstorage/out.h>

namespace flens {

template <typename FS>
std::ostream &
operator<<(std::ostream &out, const GeMatrix<FS> &A)
{
    typedef typename GeMatrix<FS>::IndexType   IndexType;
    typedef typename GeMatrix<FS>::ElementType ElementType;    
#   ifdef FLENS_IO_WITH_RANGES
    IndexType defaultIndexBase = FS::defaultIndexBase;

    out << "[";
    (A.firstRow()==defaultIndexBase) ? out << A.numRows()
                                     : out << A.firstRow()
                                           << ".."
                                           << A.lastRow();
    out << ", ";
    (A.firstCol()==defaultIndexBase) ? out << A.numCols()
                                     : out << A.firstCol()
                                           << ".."
                                           << A.lastCol();
    out << "]" << std::endl;
#   endif // FLENS_IO_WITH_RANGES

    out << std::endl;
    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            if (IsNotComplex<ElementType>::value)
                out.width(13);
            else
                out.width(28);
            out << A(i,j) << " ";
        }
        out << std::endl;
    }
    return out;
}

template <typename FS>
std::ostream &
operator<<(std::ostream &out, const HeMatrix<FS> &A)
{
    typedef typename HeMatrix<FS>::ElementType  ElementType;
    typedef typename HeMatrix<FS>::IndexType    IndexType;
#   ifdef FLENS_IO_WITH_RANGES
    IndexType defaultIndexBase = FS::defaultIndexBase;

    if ((A.firstRow()==defaultIndexBase)
     && (A.firstCol()==defaultIndexBase))
    {
        out << "[" << A.dim() << ", " << A.dim()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstRow() << ".." << A.lastRow()
            << ", " << A.firstCol() << ".." << A.lastCol()
            << "]" << std::endl;
    }
#   endif // FLENS_IO_WITH_RANGES

    out << std::endl;

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            if (IsNotComplex<ElementType>::value)
                out.width(13);
            else
                out.width(28);

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
        }
        out << std::endl;
    }
    return out;
}

template <typename FS>
std::ostream &
operator<<(std::ostream &out, const SyMatrix<FS> &A)
{
    typedef typename GeMatrix<FS>::IndexType IndexType;
    typedef typename GeMatrix<FS>::ElementType ElementType;    
#   ifdef FLENS_IO_WITH_RANGES
    IndexType defaultIndexBase = FS::defaultIndexBase;

    if ((A.firstRow()==defaultIndexBase)
     && (A.firstCol()==defaultIndexBase))
    {
        out << "[" << A.dim() << ", " << A.dim()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstRow() << ".." << A.lastRow()
            << ", " << A.firstCol() << ".." << A.lastCol()
            << "]" << std::endl;
    }
#   endif // FLENS_IO_WITH_RANGES

    out << std::endl;

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            if (IsNotComplex<ElementType>::value)
                out.width(13);
            else
                out.width(28);

            out << ((A.upLo()==cxxblas::Upper)
                    ? A(std::min(i,j), std::max(i,j))
                    : A(std::max(i,j), std::min(i,j)));
            out << " ";
        }
        out << std::endl;
    }
    return out;
}

template <typename FS>
std::ostream &
operator<<(std::ostream &out, const TrMatrix<FS> &A)
{
    typedef typename GeMatrix<FS>::IndexType    IndexType;
    typedef typename GeMatrix<FS>::ElementType  ElementType;

#   ifdef FLENS_IO_WITH_RANGES
    IndexType defaultIndexBase = FS::defaultIndexBase;

    if ((A.firstRow()==defaultIndexBase)
     && (A.firstCol()==defaultIndexBase))
    {
        out << "[" << A.dim() << ", " << A.dim()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstRow() << ".." << A.lastRow()
            << ", " << A.firstCol() << ".." << A.lastCol()
            << "]" << std::endl;
    }
#   endif // FLENS_IO_WITH_RANGES

    out << std::endl;

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            if (IsNotComplex<ElementType>::value)
                out.width(13);
            else
                out.width(28);

            if (i==j) {
                (A.diag()==cxxblas::Unit) ? out << ElementType(1)
                                          : out << A(i,j);
            } else {
                if (((i>j) && (A.upLo()==cxxblas::Lower))
                 || ((i<j) && (A.upLo()==cxxblas::Upper))) {
                    out << A(i,j);
                } else {
                    out << ElementType(0);
                }
            }
            out << " ";
        }
        out << std::endl;
    }
    return out;
}

} // namespace flens

#endif // FLENS_IO_FULLSTORAGE_OUT_TCC
