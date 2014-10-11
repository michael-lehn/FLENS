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

#ifndef FLENS_IO_FULLSTORAGE_SAVE_TCC
#define FLENS_IO_FULLSTORAGE_SAVE_TCC 1


#include <cxxstd/fstream.h>
#include <cxxstd/iomanip.h>

#include <cxxblas/typedefs.h>
#include <flens/auxiliary/iscomplex.h>
#include <flens/auxiliary/isinteger.h>
#include <flens/io/fullstorage/save.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens {

template <typename FS>
bool
save(std::string filename, const GeMatrix<FS> &A)
{

    typedef typename FS::IndexType   IndexType;
    typedef typename FS::ElementType ElementType;

    std::ofstream ofs( filename.c_str(), std::ios::binary );

    if (ofs.is_open()==false) {
        return false;
    }

    IndexType numRows = A.numRows(), numCols = A.numCols();
    IndexType firstRow = A.firstRow(), firstCol = A.firstCol();

    ofs.write(reinterpret_cast<char*>(&numRows), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&numCols), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&firstRow), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&firstCol), sizeof(IndexType));

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            ofs.write(reinterpret_cast<const char*>(&(A(i,j))),
                      sizeof(ElementType));
        }
    }

    ofs.close();
    return true;


}

template <typename FS>
bool
save(std::string filename, const HeMatrix<FS> &A)
{

    typedef typename FS::IndexType   IndexType;
    typedef typename FS::ElementType ElementType;

    std::ofstream ofs( filename.c_str(), std::ios::binary );

    if (ofs.is_open()==false) {
        return false;
    }

    IndexType dim = A.dim();
    IndexType firstIndex = A.firstRow();

    ofs.write(reinterpret_cast<char*>(&dim), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&firstIndex), sizeof(IndexType));

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=i; ++j) {
            if (A.upLo()==cxxblas::Lower) {
                ofs.write(reinterpret_cast<const char*>(&(A(i,j))),
                          sizeof(ElementType) );
            } else {
                ElementType alpha = cxxblas::conjugate(A(j,i));
                ofs.write(reinterpret_cast<const char*>(&alpha),
                          sizeof(ElementType) );
            }
        }
    }

    ofs.close();
    return true;
}

template <typename FS>
bool
save(std::string filename, const SyMatrix<FS> &A)
{

    typedef typename FS::IndexType   IndexType;
    typedef typename FS::ElementType ElementType;

    std::ofstream ofs( filename.c_str(), std::ios::binary );

    if (ofs.is_open()==false) {
        return false;
    }

    IndexType dim = A.dim();
    IndexType firstIndex = A.firstRow();

    ofs.write(reinterpret_cast<char*>(&dim), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&firstIndex), sizeof(IndexType));

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=i; ++j) {
            if (A.upLo()==cxxblas::Lower) {
                ofs.write(reinterpret_cast<const char*>(&(A(i,j))),
                          sizeof(ElementType) );
            } else {
                ofs.write(reinterpret_cast<const char*>(&(A(j,i))),
                          sizeof(ElementType) );
            }
        }
    }

    ofs.close();
    return true;
}

template <typename FS>
bool
save(std::string filename, const TrMatrix<FS> &A)
{
    typedef typename FS::IndexType   IndexType;
    typedef typename FS::ElementType ElementType;

    std::ofstream ofs( filename.c_str(), std::ios::binary );

    if (ofs.is_open()==false) {
        return false;
    }

    IndexType   numRows = A.numRows();
    IndexType   numCols = A.numCols();
    IndexType   firstRow = A.firstRow();
    IndexType   firstCol = A.firstCol();
    StorageUpLo upLo = A.upLo();
    Diag        diag = A.diag();

    ofs.write(reinterpret_cast<char*>(&numRows), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&numCols), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&firstRow), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&firstCol), sizeof(IndexType));
    ofs.write(reinterpret_cast<char*>(&upLo), sizeof(StorageUpLo));
    ofs.write(reinterpret_cast<char*>(&diag), sizeof(Diag));


    if (upLo==cxxblas::Lower) {
        for (IndexType i=A.firstRow(); i <= A.lastRow(); ++i) {
            const IndexType jmax = i-A.firstRow()+A.firstCol();
            for (IndexType j=A.firstCol(); j<jmax; ++j) {
                ofs.write(reinterpret_cast<const char*>(&(A(i,j))),
                          sizeof(ElementType) );
            }
            if (diag == cxxblas::NonUnit) {
                ofs.write(reinterpret_cast<const char*>(&(A(i,jmax))),
                          sizeof(ElementType) );
            }
        }
    } else {
        for (IndexType i=A.firstRow(); i <= A.lastRow(); ++i){
            const IndexType jmin = i-A.firstRow()+A.firstCol();
            if (diag == cxxblas::NonUnit) {
                ofs.write(reinterpret_cast<const char*>(&(A(i,jmin))),
                          sizeof(ElementType) );
            }
            for (IndexType j=jmin+1; j<=A.lastCol(); ++j) {
                ofs.write(reinterpret_cast<const char*>(&(A(i,j))),
                          sizeof(ElementType) );
            }
        }
    }

    ofs.close();
    return true;
}


template <typename FS>
bool
saveMatrixMarket(std::string filename, const GeMatrix<FS> &A,
                 std::string comment, int precision)
{
    using std::endl;
    using std::setprecision;
    using std::setw;

    typedef typename FS::IndexType   IndexType;
    typedef typename FS::ElementType ElementType;

    std::ofstream ofs( filename.c_str(), std::ios::out );

    if (ofs.is_open()==false) {
        return false;
    }

    if (IsInteger<ElementType>::value) {
        ofs << "%%MatrixMarket matrix array integer general" << endl;
    } else if (IsNotComplex<ElementType>::value) {
        ofs << "%%MatrixMarket matrix array real general" << endl;
    } else {
        ofs << "%%MatrixMarket matrix array complex general" << endl;
    }

    if (comment!="") {
        size_t j = comment.find_first_of('\n');
        while(j!=std::string::npos) {
            comment.replace(j, 1, "\n% ");
            j = comment.find_first_of('\n', j+3);
        }
        ofs << "% " << comment << endl;
    }

    ofs << A.numRows() << " " << A.numCols() << endl;

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            ofs << setw(precision+6) << setprecision(precision)
                << cxxblas::real(A(i,j));
            if (IsComplex<ElementType>::value) {
                 ofs << setw(precision+7) << setprecision(precision)
                     << cxxblas::imag(A(i,j));
            }
            ofs << endl;;
        }
    }
    ofs.close();
    return true;

}

template <typename FS>
bool
saveMatrixMarket(std::string filename, const SyMatrix<FS> &A,
                 std::string comment, int precision)
{
    using std::endl;
    using std::setprecision;
    using std::setw;

    typedef typename FS::IndexType   IndexType;
    typedef typename FS::ElementType ElementType;

    std::ofstream ofs(filename.c_str(), std::ios::out);

    if (ofs.is_open()==false)
        return false;

    if (IsInteger<ElementType>::value) {
        ofs << "%%MatrixMarket matrix array integer symmetric" << endl;
    } else if (IsNotComplex<ElementType>::value) {
        ofs << "%%MatrixMarket matrix array real symmetric" << endl;
    } else {
        ofs << "%%MatrixMarket matrix array complex symmetric" << endl;
    }

    if (comment!="") {
        size_t j = comment.find_first_of('\n');
        while (j!=std::string::npos) {
            comment.replace(j, 1, "\n% ");
            j = comment.find_first_of('\n', j+3);
        }
        ofs << "% " << comment << endl;
    }

    ofs << A.numRows() << " " << A.numCols() << endl;

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j = i; j <= A.lastCol(); ++j) {
            if (A.upLo()==Upper) {
                ofs << setw(precision+6) << setprecision(precision)
                    << cxxblas::real(A(i,j));
                if (IsComplex<ElementType>::value) {
                    ofs << setw(precision+7) << setprecision(precision)
                        << cxxblas::imag(A(i,j));
                }
            } else {
                ofs << setw(precision+6) << setprecision(precision)
                    <<  cxxblas::real(A(j,i));
                if (IsComplex<ElementType>::value) {
                    ofs << setw(precision+7) << setprecision(precision)
                        << cxxblas::imag(A(j,i));
                }
            }
            ofs << endl;;
        }
    }
    ofs.close();
    return true;
}

template <typename FS>
typename RestrictTo<IsComplex<typename FS::ElementType>::value, bool>::Type
saveMatrixMarket(std::string filename, const HeMatrix<FS> &A,
                 std::string comment, int precision)
{
    using std::endl;
    using std::setprecision;
    using std::setw;

    typedef typename FS::IndexType   IndexType;

    std::ofstream ofs(filename.c_str(), std::ios::out);

    if (ofs.is_open()==false) {
        return false;
    }

    ofs << "%%MatrixMarket matrix array complex hermitian" << endl;

    if (comment!="") {
        size_t j = comment.find_first_of('\n');
        while(j!=std::string::npos) {
            comment.replace(j, 1, "\n% ");
            j = comment.find_first_of('\n', j+3);
        }
        ofs << "% " << comment << endl;
    }

    ofs << A.numRows() << " " << A.numCols() << endl;

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=i; j<=A.lastCol(); ++j) {
            if (A.upLo()==Upper) {
                ofs << setw(precision+6) << setprecision(precision)
                    <<  cxxblas::real(A(i,j));
                if (i==j) {
                    ofs << setw(precision+7) << setprecision(precision) << 0;
                } else {
                    ofs << setw(precision+7) << setprecision(precision)
                        << -cxxblas::imag(A(i,j));
                }

            } else {
                ofs << setw(precision+6) << setprecision(precision)
                    <<  cxxblas::real(A(j,i));
                if (i==j) {
                    ofs << setw(precision+7) << setprecision(precision) << 0;
                } else {
                    ofs << setw(precision+7) << setprecision(precision)
                        << cxxblas::imag(A(j,i));
                }
            }
            ofs << endl;;
        }
    }
    ofs.close();
    return true;
}

} // namespace flens

#endif // FLENS_IO_FULLSTORAGE_SAVE_TCC

