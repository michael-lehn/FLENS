/*
 *   Copyright (c) 2007, Michael Lehn
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

#ifndef FLENS_PACKEDSTORAGE_H
#define FLENS_PACKEDSTORAGE_H 1

#include <flens/array.h>

namespace flens {

// forward decalarations
template <typename T, StorageOrder Order, StorageUpLo UpLo>
    class PackedStorage;

template <typename T, StorageOrder Order, StorageUpLo UpLo>
    class PackedStorageView;

template <typename T, StorageOrder Order, StorageUpLo UpLo>
class ConstPackedStorageView
{
    public:
        typedef T                                        ElementType;
        typedef ConstPackedStorageView<T, Order, UpLo>   ConstView;
        typedef PackedStorageView<T, Order, UpLo>        View;
        typedef PackedStorage<T, Order, UpLo>            NoView;

        ConstPackedStorageView(const ConstArrayView<T> &data,
                               int dim, int indexBase=1);

        ConstPackedStorageView(const ArrayView<T> &data,
                               int dim, int indexBase=1);

        ConstPackedStorageView(const ConstPackedStorageView<T,Order,UpLo> &rhs);

        ~ConstPackedStorageView();

        //-- operators ---------------------------------------------------------

        const T &
        operator()(int row, int col) const;

        //-- methods -----------------------------------------------------------

        int
        firstIndex() const;

        int
        lastIndex() const;

        int
        dim() const;

        const T *
        data() const;

        ConstView
        view() const;

    private:
        ConstArrayView<T>  _data;
        int                _dim, _indexBase;
};

//------------------------------------------------------------------------------

template <typename T, StorageOrder Order, StorageUpLo UpLo>
struct StorageInfo<ConstPackedStorageView<T, Order, UpLo> >
{
    static const StorageOrder order = Order;
    static const StorageUpLo upLo = UpLo;
};

//------------------------------------------------------------------------------

template <typename T, StorageOrder Order, StorageUpLo UpLo>
class PackedStorageView
{
    public:
        typedef T                                        ElementType;
        typedef ConstPackedStorageView<T, Order, UpLo>   ConstView;
        typedef PackedStorageView<T, Order, UpLo>        View;
        typedef PackedStorage<T, Order, UpLo>            NoView;

        PackedStorageView(const ArrayView<T> &data,
                          int dim, int indexBase=1);

        PackedStorageView(const PackedStorageView<T, Order, UpLo> &rhs);

        ~PackedStorageView();

        //-- operators ---------------------------------------------------------

        PackedStorageView<T, Order, UpLo> &
        operator=(const NoView &rhs);

        PackedStorageView<T, Order, UpLo> &
        operator=(const View &rhs);

        PackedStorageView<T, Order, UpLo> &
        operator=(const ConstView &rhs);

        const T &
        operator()(int row, int col) const;

        T &
        operator()(int row, int col);

        //-- methods -----------------------------------------------------------

        int
        firstIndex() const;

        int
        lastIndex() const;

        int
        dim() const;

        const T *
        data() const;

        T *
        data();

        ConstView
        view() const;

        View
        view();

    private:
        ArrayView<T>  _data;
        int           _dim, _indexBase;
};

//------------------------------------------------------------------------------

template <typename T, StorageOrder Order, StorageUpLo UpLo>
struct StorageInfo<PackedStorageView<T, Order, UpLo> >
{
    static const StorageOrder order = Order;
    static const StorageUpLo upLo = UpLo;
};

//------------------------------------------------------------------------------

template <typename T, StorageOrder Order, StorageUpLo UpLo>
class PackedStorage
{
    public:
        typedef T                                        ElementType;
        typedef ConstPackedStorageView<T, Order, UpLo>   ConstView;
        typedef PackedStorageView<T, Order, UpLo>        View;
        typedef PackedStorage<T, Order, UpLo>            NoView;

        PackedStorage();

        PackedStorage(int dim, int indexBase=1);

        PackedStorage(const PackedStorage<T, Order, UpLo> &rhs);

        ~PackedStorage();

        //-- operators ---------------------------------------------------------

        PackedStorage<T, Order, UpLo> &
        operator=(const NoView &rhs);

        PackedStorage<T, Order, UpLo> &
        operator=(const View &rhs);

        PackedStorage<T, Order, UpLo> &
        operator=(const ConstView &rhs);

        const T &
        operator()(int row, int col) const;

        T &
        operator()(int row, int col);

        //-- methods -----------------------------------------------------------

        int
        firstIndex() const;

        int
        lastIndex() const;

        int
        dim() const;

        const T *
        data() const;

        T *
        data();

        void
        resize(int dim, int indexBase=1);

        ConstView
        view() const;

        View
        view();

    private:
        Array<T>  _data;
        int       _dim, _indexBase;
};

//------------------------------------------------------------------------------

template <typename T, StorageOrder Order, StorageUpLo UpLo>
struct StorageInfo<PackedStorage<T, Order, UpLo> >
{
    static const StorageOrder order = Order;
    static const StorageUpLo upLo = UpLo;
};


} // namespace flens

#include <flens/packedstorage.tcc>

#endif // FLENS_PACKEDSTORAGE_H
