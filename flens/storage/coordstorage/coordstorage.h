/*
 *   Copyright (c) 2007-2012, Michael Lehn
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

#ifndef FLENS_STORAGE_COORDSTORAGE_COORDSTORAGE_H
#define FLENS_STORAGE_COORDSTORAGE_COORDSTORAGE_H 1

#include <cxxstd/algorithm.h>
#include <cxxstd/cstddef.h>
#include <cxxstd/vector.h>
#include <flens/storage/indexoptions.h>

namespace flens {

// forward declarations
template <typename T, typename IndexType>
    struct Coord;

template <typename T>
    struct CoordProxy;

struct CoordRowColCmp;


template <typename T,
          typename Cmp = CoordRowColCmp,
          typename I = IndexOptions<> >
class CoordStorage
{
    public:
        typedef T                        ElementType;
        typedef typename I::IndexType    IndexType;
        typedef CoordProxy<T>            ElementProxy;

        typedef Coord<T, IndexType>      CoordType;
        typedef std::vector<CoordType>   CoordVector;

        static const IndexType                        defaultIndexBase
                                                          = I::defaultIndexBase;

        CoordStorage(IndexType numRows, IndexType numCols,
                     IndexType densityEstimate = 1,
                     IndexType indexBase = I::defaultIndexBase);

        ~CoordStorage();

        //-- operators ---------------------------------------------------------

        ElementProxy
        operator()(IndexType row, IndexType col);

        //-- methods -----------------------------------------------------------

        const IndexType
        indexBase() const;

        const IndexType
        firstRow() const;

        const IndexType
        lastRow() const;

        const IndexType
        firstCol() const;

        const IndexType
        lastCol() const;

        const IndexType
        numRows() const;

        const IndexType
        numCols() const;

        void
        accumulate() const;

        const CoordVector &
        coordVector() const;

        const IndexType
        numNonZeros() const;

    //private:
        // Lehn:  I do not allow copying matrices with coordiante storage unless
        //        someone gives me a reason.
        CoordStorage(const CoordStorage &rhs);

        IndexType  numRows_, numCols_;
        IndexType  indexBase_;

        mutable CoordVector  coord_;
        mutable size_t       lastSortedCoord_;
        mutable bool         isSorted_;
        mutable bool         isAccumulated_;
        Cmp                  less_;
};

//-- Coord ---------------------------------------------------------------------

template <typename T, typename IndexType>
struct Coord
{
    Coord(IndexType row, IndexType col, const T &value);

    IndexType   row, col;
    T           value;
};

//-- CoordProxy ----------------------------------------------------------------

template <typename T>
struct CoordProxy
{
    CoordProxy(T *value);

    void
    operator+=(const T &x);

    void
    operator-=(const T &x);

    private:
        T  *value;

//
//      We only allow accumulation of data.  Assignment would be too expensive.
//
        void
        operator=(const T &x);
};

//-- CoordRowColCmp ------------------------------------------------------------
//
//  This ordering gets used for converting the coordinate storage format to
//  the compressed row storage format
//
struct CoordRowColCmp
{
        // return true if a < b
        //   <=>  a.row<b.row  or  a.row==b.row and a.col<b.col
        template <typename T, typename IndexType>
            bool
            operator()(const Coord<T, IndexType> &a,
                       const Coord<T, IndexType> &b) const;
};

//-- CoordColRowCmp ------------------------------------------------------------
//
//  This ordering gets used for converting the coordinate storage format to
//  the compressed row storage format
//
struct CoordColRowCmp
{
        // return true if a < b
        //   <=>  a.col<b.col  or  a(2)==b(2) and a.row<b.row
        template <typename T, typename IndexType>
            bool
            operator()(const Coord<T, IndexType> &a,
                       const Coord<T, IndexType> &b) const;
};

} // namespace flens

#endif // FLENS_STORAGE_COORDSTORAGE_COORDSTORAGE_H
