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

#ifndef PLAYGROUND_FLENS_MPI_REDUCE_REDUCE_TCC
#define PLAYGROUND_FLENS_MPI_REDUCE_REDUCE_TCC 1

#include<playground/flens/mpi/mpi-flens.h>

namespace flens { namespace mpi {

#ifdef WITH_MPI

template <typename T>
typename RestrictTo<IsReal<T>::value ||
                    IsInteger<T>::value,
                    T>::Type
MPI_reduce_max(const T &x, const int root, const MPI::Comm &communicator)
{
    using namespace MPI;

    typedef typename MPI_Type<T>::PrimitiveType  PT;

    T sum(0);
    communicator.Reduce(reinterpret_cast<const PT *>(&x),
                        reinterpret_cast<PT *>(&sum),
                        MPI_Type<T>::size, MPI_Type<T>::Type(), MPI_MAX, 0);
    return sum;

}

template <typename T>
typename RestrictTo<IsReal<T>::value ||
                    IsInteger<T>::value,
                    T>::Type
MPI_reduce_min(const T &x, const int root, const MPI::Comm &communicator)
{
    using namespace MPI;

    typedef typename MPI_Type<T>::PrimitiveType  PT;

    T sum(0);
    communicator.Reduce(reinterpret_cast<const PT *>(&x),
                        reinterpret_cast<PT *>(&sum),
                        MPI_Type<T>::size, MPI_Type<T>::Type(), MPI_MIN, 0);
    return sum;

}


template <typename T>
typename RestrictTo<MPI_Type<T>::Compatible,
                    T>::Type
MPI_reduce_sum(const T &x, const int root, const MPI::Comm &communicator)
{
    using namespace MPI;

    typedef typename MPI_Type<T>::PrimitiveType  PT;

    T sum(0);
    communicator.Reduce(reinterpret_cast<const PT *>(&x),
                        reinterpret_cast<PT *>(&sum),
                        MPI_Type<T>::size, MPI_Type<T>::Type(), MPI_SUM, 0);
    return sum;

}

template <typename VX, typename VSUM>
typename RestrictTo<IsDenseVector<VX>::value &&
                    IsDenseVector<VSUM>::value,
                    void>::Type
MPI_reduce_sum(VX &&x, VSUM &&sum, const int root,
               const MPI::Comm &communicator)
{
    using namespace MPI;

    typedef typename MPI_Type<T>::PrimitiveType  PT;

    typedef typename RemoveRef<VX>::Type   VectorX;
    typedef typename VectorX::ElementType T;
    typedef typename VectorX::IndexType   IndexType;

    const int rank = MPI_rank();

    ASSERT( sum.length()==x.length() || root!=MPI_rank() );
    ASSERT( root<MPI_size() );

#ifndef NDEBUG
    // Check corrent length and stride
    IndexType length = x.length();
    IndexType stride = x.stride();
    communicator.Bcast(&length, MPI_Type<IndexType>::size,
                       MPI_Type<IndexType>::Type(), root);
    communicator.Bcast(&stride, MPI_Type<IndexType>::size,
                       MPI_Type<IndexType>::Type(), root);

    ASSERT( x.length()==length );
    ASSERT( x.stride()==stride );
#endif

    if ( x.stride()==1
     && ((rank==root && sum.stride()==x.stride()) || rank!=root ))
    {

        T *psum = NULL;
        if (root == rank) {
            psum = sum.data();
        }

        communicator.Reduce(reinterpret_cast<PT *>(x.data()),
                            reinterpret_cast<PT *>(psum),
                            x.length()*MPI_Type<T>::size,
                            MPI_Type<T>::Type(), MPI_SUM, 0);

    } else {

        const IndexType i0 = x.firstIndex(),
                        i1 = x.lastIndex();
        const IndexType j0 = sum.firstIndex();

        for (IndexType i=i0, j=j0; i<=i1; ++i, ++j) {

            T *psum = NULL;
            if (root == rank) {
                psum = &sum(j);
            }

            communicator.Reduce(reinterpret_cast<PT *>(&x(i)),
                                reinterpret_cast<PT *>(psum),
                                MPI_Type<T>::size, MPI_Type<T>::Type(),
                                MPI_SUM, 0);


        }
    }

}


template <typename MA, typename MSUM>
typename RestrictTo<IsGeMatrix<MA>::value &&
                    IsGeMatrix<MSUM>::value,
                    void>::Type
MPI_reduce_sum(MA &&A, MSUM &&Sum, const int root,
               const MPI::Comm &communicator)
{

    using namespace MPI;

    typedef typename RemoveRef<MA>::Type   MatrixA;
    typedef typename MatrixA::ElementType T;
    typedef typename MatrixA::IndexType   IndexType;

    const Underscore<IndexType> _;

    const int rank = MPI_rank();

    if ( Sum.numRows()==0 && Sum.numCols()==0 && root==MPI_rank()) {
        Sum.resize(A.numRows(), A.numCols());
    }

    ASSERT( (Sum.numRows()==A.numRows() &&
            (Sum.numCols()==A.numCols()) &&
            (A.order()==Sum.order() ))
          || root!=MPI_rank() );
    ASSERT( root<MPI_size() );

#ifndef NDEBUG
    // Check correct size
    IndexType numCols = A.numCols();
    IndexType numRows = A.numRows();

    communicator.Bcast(&numRows, MPI_Type<IndexType>::size,
                       MPI_Type<IndexType>::Type(), root);
    communicator.Bcast(&numCols, MPI_Type<IndexType>::size,
                       MPI_Type<IndexType>::Type(), root);

    ASSERT( A.numRows()==numRows );
    ASSERT( A.numCols()==numCols );
#endif

    if ( A.order()==ColMajor ) {
        const IndexType i0 = A.firstCol(),
                        i1 = A.lastCol();
        const IndexType j0 = Sum.firstCol();

        for (IndexType i=i0, j=j0; i<=i1; ++i, ++j) {
            auto x   = A(_, i);
            auto sum = Sum(_, j);
            MPI_reduce_sum(x, sum, root, communicator);
        }
    } else
        const IndexType i0 = A.firstRow(),
                        i1 = A.lastRow();
        const IndexType j0 = Sum.firstRow();

        for (IndexType i=i0, j=j0; i<=i1; ++i, ++j) {
            auto x   = A(i,_);
            auto sum = Sum(j,_);
            MPI_reduce_sum(x, sum, root, communicator);
        }
    }

}

#else

//--- Max ---------------------------------------------------------------------
template <typename T>
typename RestrictTo<IsReal<T>::value ||
                    IsInteger<T>::value,
                    T>::Type
MPI_reduce_max(const T &x, const int root)
{
    return x;
}

//--- Min ---------------------------------------------------------------------
template <typename T>
typename RestrictTo<IsReal<T>::value ||
                    IsInteger<T>::value,
                    T>::Type
MPI_reduce_min(const T &x, const int root)
{
    return x;
}

//--- Sum ---------------------------------------------------------------------
template <typename T>
T
MPI_reduce_sum(const T &x, const int root)
{
    ASSERT( root==0 );
    return x;
}

template <typename X, typename SUM>
typename RestrictTo<(IsDenseVector<X>::value &&
                     IsDenseVector<SUM>::value) ||
                    (IsGeMatrix<X>::value &&
                     IsGeMatrix<SUM>::value),
                     void>::Type
MPI_reduce_sum(X &&x, SUM &&sum, const int root)
{
    ASSERT( root==0 );
    sum = x;
}

#endif // WITH_MPI

} }

#endif // PLAYGROUND_FLENS_MPI_REDUCE_REDUCE_TCC
