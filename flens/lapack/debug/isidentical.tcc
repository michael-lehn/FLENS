/*
 *   Copyright (c) 2011, Michael Lehn
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

#ifndef FLENS_LAPACK_DEBUG_ISIDENTICAL_TCC
#define FLENS_LAPACK_DEBUG_ISIDENTICAL_TCC 1

#include <complex>
#include <flens/lapack/debug/hex.h>
#include <flens/lapack/debug/isidentical.h>

namespace flens { namespace lapack {

template <typename X, typename Y>
bool
isDifferent(const X &x, const Y &y)
{
    using std::isnan;

    if (isnan(x) && isnan(y)) {
        return false;
    }
    return x!=y;
}

template <typename X, typename Y>
bool
isDifferent(const std::complex<X> &x, const std::complex<Y> &y)
{
    using std::isnan;

    return isDifferent(x.real(), y.real()) || isDifferent(x.imag(), y.imag());
}

template <typename X, typename Y>
bool
isIdentical(const X &x, const Y &y, const char *xName, const char *yName)
{
    if (isDifferent(x,y)) {
        std::cerr.precision(50);
        std::cerr << xName << " = " << x
                  << std::endl
                  << yName << " = " << y
                  << std::endl
                  << xName << " - " << yName << " = "
                  << x - y
                  << std::endl;
        std::cerr << "hex(" << xName << ") = " << hex(x)
                  << std::endl
                  << "hex(" << yName << ") = " << hex(y)
                  << std::endl
                  << "hex(" << xName << " - " << yName << ") = "
                  << hex(x - y)
                  << std::endl;
        return false;
    }
    return true;
}

template <typename VX, typename VY>
bool
isIdentical(const DenseVector<VX> &x, const DenseVector<VY> &y,
            const char *xName, const char *yName)
{
    typedef typename DenseVector<VX>::IndexType IndexType;

    if (x.length()!=y.length()) {
        std::cerr << xName << ".length() = " << x.length() << ", "
                  << yName << ".length() = " << y.length()
                  << std::endl;
        return false;
    }
    if (x.firstIndex()!=y.firstIndex()) {
        std::cerr << xName << ".firstIndex() = " << x.firstIndex() << ", "
                  << yName << ".firstIndex() = " << y.firstIndex()
                  << std::endl;
        return false;
    }
    if (x.endIndex()!=y.endIndex()) {
        std::cerr << xName << ".endIndex() = " << x.endIndex() << ", "
                  << yName << ".endIndex() = " << y.endIndex()
                  << std::endl;
        return false;
    }
    if (x.inc()!=y.inc()) {
        std::cerr << xName << ".inc() = " << x.inc() << ", "
                  << yName << ".inc() = " << y.inc()
                  << std::endl;
        return false;
    }

    for (IndexType i=x.firstIndex(); i!=x.endIndex(); i+=x.inc()) {
        if (isDifferent(x(i), y(i))) {
            std::cerr.precision(50);
            std::cerr << xName << "(" << i << ") = " << x(i)
                      << std::endl
                      << yName << "(" << i << ") = " << y(i)
                      << std::endl
                      << xName << "(" << i << ") - "
                      << yName << "(" << i << ") = "
                      << x(i) - y(i)
                      << std::endl;
            std::cerr << "hex(" << xName << "(" << i << " )) = " << hex(x(i))
                      << std::endl
                      << "hex(" << yName << "(" << i << " )) = " << hex(y(i))
                      << std::endl
                      << "hex(" << xName << "(" << i << ") - "
                                << yName << "(" << i << ")) = "
                      << hex(x(i) - y(i))
                      << std::endl;
            return false;
        }
    }
    return true;
}

template <typename MA, typename MB>
bool
isIdentical(const GeMatrix<MA> &A, const GeMatrix<MB> &B,
            const char *AName, const char *BName)
{
    typedef typename GeMatrix<MA>::IndexType IndexType;

    if (A.numRows()*A.numCols()==0 && B.numRows()*B.numCols()==0) {
        return true;
    }

    if (A.numRows()!=B.numRows()) {
        std::cerr << AName << ".numRows() = " << A.numRows() << ", "
                  << BName << ".numRows() = " << B.numRows()
                  << std::endl;
        return false;
    }
    if (A.numCols()!=B.numCols()) {
        std::cerr << AName << ".numCols() = " << A.numCols() << ", "
                  << BName << ".numCols() = " << B.numCols()
                  << std::endl;
        return false;
    }
    if (A.firstRow()!=B.firstRow()) {
        std::cerr << AName << ".firstRow() = " << A.firstRow() << ", "
                  << BName << ".firstRow() = " << B.firstRow()
                  << std::endl;
        return false;
    }
    if (A.firstCol()!=B.firstCol()) {
        std::cerr << AName << ".firstCol() = " << A.firstCol() << ", "
                  << BName << ".firstCol() = " << B.firstCol()
                  << std::endl;
        return false;
    }
    if (A.lastRow()!=B.lastRow()) {
        std::cerr << AName << ".lastRow() = " << A.lastRow() << ", "
                  << BName << ".lastRow() = " << B.lastRow()
                  << std::endl;
        return false;
    }
    if (A.lastCol()!=B.lastCol()) {
        std::cerr << AName << ".lastCol() = " << A.lastCol() << ", "
                  << BName << ".lastCol() = " << B.lastCol()
                  << std::endl;
        return false;
    }

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            if (isDifferent(A(i,j), B(i,j))) {
                std::cerr.precision(50);
                std::cerr << AName << "(" << i << ", " << j << ") = " << A(i,j)
                          << std::endl
                          << BName << "(" << i << ", " << j << ") = " << B(i,j)
                          << std::endl
                          << AName << "(" << i << ", " << j << ") - "
                          << BName << "(" << i << ", " << j << ") = "
                          << A(i,j) - B(i,j) << std::endl;
                std::cerr << "hex(" << AName << "(" << i << ", " << j << ")) = "
                          << hex(A(i,j)) << std::endl
                          << "hex(" << BName << "(" << i << ", " << j << ")) = " 
                          << hex(B(i,j)) << std::endl
                          << "hex(" << AName << "(" << i << ", " << j << ")) - "
                          << "hex(" << BName << "(" << i << ", " << j << ")) = "
                          << hex(A(i,j) - B(i,j)) << std::endl;
                 return false;
            }
        }
    }
    return true;
}

template <typename MA, typename MB>
bool
isIdentical(const HeMatrix<MA> &A, const HeMatrix<MB> &B,
            const char *AName, const char *BName)
{
    typedef typename SyMatrix<MA>::IndexType IndexType;

    if (A.dim()==0 && B.dim()==0) {
        return true;
    }

    if (A.dim()!=B.dim()) {
        std::cerr << AName << ".dim() = " << A.dim() << ", "
                  << BName << ".dim() = " << B.dim()
                  << std::endl;
        return false;
    }
    if (A.firstRow()!=B.firstRow()) {
        std::cerr << AName << ".firstRow() = " << A.firstRow() << ", "
                  << BName << ".firstRow() = " << B.firstRow()
                  << std::endl;
        return false;
    }
    if (A.firstCol()!=B.firstCol()) {
        std::cerr << AName << ".firstCol() = " << A.firstCol() << ", "
                  << BName << ".firstCol() = " << B.firstCol()
                  << std::endl;
        return false;
    }
    if (A.lastRow()!=B.lastRow()) {
        std::cerr << AName << ".lastRow() = " << A.lastRow() << ", "
                  << BName << ".lastRow() = " << B.lastRow()
                  << std::endl;
        return false;
    }
    if (A.lastCol()!=B.lastCol()) {
        std::cerr << AName << ".lastCol() = " << A.lastCol() << ", "
                  << BName << ".lastCol() = " << B.lastCol()
                  << std::endl;
        return false;
    }

    if (A.upLo()!=B.upLo()) {
        std::cerr << AName << ".upLo() = " << A.upLo() << ", "
                  << BName << ".upLo() = " << B.upLo()
                  << std::endl;
        return false;
    }

    bool failed = false;
    bool isUpper = (A.upLo()==Upper);

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            if (isUpper && j>=i) {
                failed = isDifferent(A(i,j), B(i,j));
            }
            if (!isUpper && j<=i) {
                failed = isDifferent(A(i,j), B(i,j));
            }
            if (failed) {
                std::cerr.precision(50);
                std::cerr << AName << "(" << i << ", " << j << ") = "
                          << A(i,j) << std::endl
                          << BName << "(" << i << ", " << j << ") = "
                          << B(i,j) << std::endl
                          << AName << "(" << i << ", " << j << ") - "
                          << BName << "(" << i << ", " << j << ") = "
                          << A(i,j) - B(i,j)
                          << std::endl;
                std::cerr << "hex(" << AName << "(" << i << ", " << j << ")) = "
                          << hex(A(i,j)) << std::endl
                          << "hex(" << BName << "(" << i << ", " << j << ")) = "
                          << hex(B(i,j)) << std::endl
                          << "hex(" << AName << "(" << i << ", " << j << ")) - "
                          << "hex(" << BName << "(" << i << ", " << j << ")) = "
                          << hex(A(i,j) - B(i,j))
                          << std::endl;
                 return false;
            }
        }
    }
    return true;
}



template <typename MA, typename MB>
bool
isIdentical(const TrMatrix<MA> &A, const TrMatrix<MB> &B,
            const char *AName, const char *BName)
{
    typedef typename TrMatrix<MA>::IndexType IndexType;

    if (A.dim()==0 && B.dim()==0) {
        return true;
    }

    if (A.dim()!=B.dim()) {
        std::cerr << AName << ".dim() = " << A.dim() << ", "
                  << BName << ".dim() = " << B.dim()
                  << std::endl;
        return false;
    }
    if (A.firstRow()!=B.firstRow()) {
        std::cerr << AName << ".firstRow() = " << A.firstRow() << ", "
                  << BName << ".firstRow() = " << B.firstRow()
                  << std::endl;
        return false;
    }
    if (A.firstCol()!=B.firstCol()) {
        std::cerr << AName << ".firstCol() = " << A.firstCol() << ", "
                  << BName << ".firstCol() = " << B.firstCol()
                  << std::endl;
        return false;
    }
    if (A.lastRow()!=B.lastRow()) {
        std::cerr << AName << ".lastRow() = " << A.lastRow() << ", "
                  << BName << ".lastRow() = " << B.lastRow()
                  << std::endl;
        return false;
    }
    if (A.lastCol()!=B.lastCol()) {
        std::cerr << AName << ".lastCol() = " << A.lastCol() << ", "
                  << BName << ".lastCol() = " << B.lastCol()
                  << std::endl;
        return false;
    }

    if (A.upLo()!=B.upLo()) {
        std::cerr << AName << ".upLo() = " << A.upLo() << ", "
                  << BName << ".upLo() = " << B.upLo()
                  << std::endl;
        return false;
    }

    bool failed = false;
    bool isUpper = (A.upLo()==Upper);
    bool isUnit = (A.diag()==Unit);

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            if (isUnit && i==j) {
                continue;
            }
            if (isUpper && j>=i) {
                failed = isDifferent(A(i,j), B(i,j));
            }
            if (!isUpper && j<=i) {
                failed = isDifferent(A(i,j), B(i,j));
            }
            if (failed) {
                std::cerr.precision(50);
                std::cerr << AName << "(" << i << ", " << j << ") = "
                          << A(i,j) << std::endl
                          << BName << "(" << i << ", " << j << ") = "
                          << B(i,j) << std::endl
                          << AName << "(" << i << ", " << j << ") - "
                          << BName << "(" << i << ", " << j << ") = "
                          << A(i,j) - B(i,j)
                          << std::endl;
                std::cerr << "hex(" << AName << "(" << i << ", " << j << ")) = "
                          << hex(A(i,j)) << std::endl
                          << "hex(" << BName << "(" << i << ", " << j << ")) = "
                          << hex(B(i,j)) << std::endl
                          << "hex(" << AName << "(" << i << ", " << j << ")) - "
                          << "hex(" << BName << "(" << i << ", " << j << ")) = "
                          << hex(A(i,j) - B(i,j))
                          << std::endl;
                return false;
            }
        }
    }
    return true;
}

template <typename MA, typename MB>
bool
isIdentical(const SyMatrix<MA> &A, const SyMatrix<MB> &B,
            const char *AName, const char *BName)
{
    typedef typename SyMatrix<MA>::IndexType IndexType;

    if (A.dim()==0 && B.dim()==0) {
        return true;
    }

    if (A.dim()!=B.dim()) {
        std::cerr << AName << ".dim() = " << A.dim() << ", "
                  << BName << ".dim() = " << B.dim()
                  << std::endl;
        return false;
    }
    if (A.firstRow()!=B.firstRow()) {
        std::cerr << AName << ".firstRow() = " << A.firstRow() << ", "
                  << BName << ".firstRow() = " << B.firstRow()
                  << std::endl;
        return false;
    }
    if (A.firstCol()!=B.firstCol()) {
        std::cerr << AName << ".firstCol() = " << A.firstCol() << ", "
                  << BName << ".firstCol() = " << B.firstCol()
                  << std::endl;
        return false;
    }
    if (A.lastRow()!=B.lastRow()) {
        std::cerr << AName << ".lastRow() = " << A.lastRow() << ", "
                  << BName << ".lastRow() = " << B.lastRow()
                  << std::endl;
        return false;
    }
    if (A.lastCol()!=B.lastCol()) {
        std::cerr << AName << ".lastCol() = " << A.lastCol() << ", "
                  << BName << ".lastCol() = " << B.lastCol()
                  << std::endl;
        return false;
    }

    if (A.upLo()!=B.upLo()) {
        std::cerr << AName << ".upLo() = " << A.upLo() << ", "
                  << BName << ".upLo() = " << B.upLo()
                  << std::endl;
        return false;
    }

    bool failed = false;
    bool isUpper = (A.upLo()==Upper);

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            if (isUpper && j>=i) {
                failed = isDifferent(A(i,j), B(i,j));
            }
            if (!isUpper && j<=i) {
                failed = isDifferent(A(i,j), B(i,j));
            }
            if (failed) {
                std::cerr.precision(50);
                std::cerr << AName << "(" << i << ", " << j << ") = "
                          << A(i,j) << std::endl
                          << BName << "(" << i << ", " << j << ") = "
                          << B(i,j) << std::endl
                          << AName << "(" << i << ", " << j << ") - "
                          << BName << "(" << i << ", " << j << ") = "
                          << A(i,j) - B(i,j)
                          << std::endl;
                std::cerr << "hex(" << AName << "(" << i << ", " << j << ")) = "
                          << hex(A(i,j)) << std::endl
                          << "hex(" << BName << "(" << i << ", " << j << ")) = "
                          << hex(B(i,j)) << std::endl
                          << "hex(" << AName << "(" << i << ", " << j << ")) - "
                          << "hex(" << BName << "(" << i << ", " << j << ")) = "
                          << hex(A(i,j) - B(i,j))
                          << std::endl;
                 return false;
            }
        }
    }
    return true;
}


} } // namespace lapack, flens

#endif // FLENS_LAPACK_DEBUG_ISIDENTICAL_TCC
