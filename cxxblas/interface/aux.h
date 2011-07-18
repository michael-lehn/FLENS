#ifndef CXXBLAS_INTERFACE_AUX_H
#define CXXBLAS_INTERFACE_AUX_H 1

#include <cxxblas/cxxblas.cxx>

using cxxblas::StorageOrder;
using cxxblas::ColMajor;
using cxxblas::RowMajor;
using cxxblas::StorageUpLo;
using cxxblas::Upper;
using cxxblas::Lower;


template <typename IndexType, typename MA>
void
switchFullStorageOrder(StorageOrder order, IndexType m, IndexType n,
                       const MA *&A, IndexType &ldA)
{
    assert(order==ColMajor);

    IndexType _ldA = n+1;
    MA *_A = new MA[_ldA*m];

    for (IndexType i=0; i<m; ++i) {
        for (IndexType j=0; j<n; ++j) {
            _A[i*_ldA+j] = A[j*ldA+i];
        }
    }

    A = _A;
    ldA = _ldA;
}

template <typename IndexType, typename MA>
void
switchFullStorageOrder(StorageOrder orderA, IndexType m, IndexType n,
                       const MA *A, IndexType &ldA,
                       MA *&_A, IndexType &_ldA)
{
    if (orderA==RowMajor) {
        for (IndexType i=0; i<m; ++i) {
            for (IndexType j=0; j<n; ++j) {
                _A[j*_ldA+i] = A[i*ldA+j];
            }
        }
    } else {
        for (IndexType i=0; i<m; ++i) {
            for (IndexType j=0; j<n; ++j) {
                _A[i*_ldA+j] = A[j*ldA+i];
            }
        }
    }
}

template <typename IndexType, typename MA>
void
switchBandStorageOrder(StorageOrder order, IndexType m, IndexType n,
                       IndexType kl, IndexType ku,
                       const MA *&A, IndexType &ldA)
{
    assert(order==ColMajor);

    IndexType _ldA = kl+ku+1;
    MA *_A = new MA[(kl+ku+1)*m];

    for (IndexType i=0; i<m; ++i) {
        for (IndexType j=0; j<n; ++j) {
            if ((j-i<=ku) && (i-j<=kl)) {
                _A[(kl+j-i)+_ldA*i] = A[(ku+i-j)+ldA*j];
            }
        }
    }

    A = _A;
    ldA = _ldA;
}

template <typename IndexType, typename MA>
void
switchBandStorageOrder(StorageOrder order, StorageUpLo upLo,
                       IndexType n, IndexType k,
                       const MA *&A, IndexType &ldA)
{
    IndexType kl = (upLo==Lower) ? k : 0;
    IndexType ku = (upLo==Upper) ? k : 0;
    switchBandStorageOrder(order, n, n, kl, ku, A, ldA);
}

template <typename IndexType, typename MA>
void
switchPackedStorageOrder(StorageOrder order, StorageUpLo &upLo,
                         IndexType n, const MA *&A)
{
    assert(order==ColMajor);
    
    MA *_A = new MA[n*(n+1)/2];
    if (upLo==Upper) {
        for (IndexType i=0; i<n; ++i) {
            for (IndexType j=i; j<n; ++j) {
                _A[j+i*(2*n-i-1)/2] = A[i+j*(j+1)/2];
            }
        }
    } else {
        for (IndexType i=0; i<n; ++i) {
            for (IndexType j=0; j<=i; ++j) {
                _A[j+i*(i+1)/2] = A[i+j*(2*n-j-1)/2];
            }
        }
    }
    A = _A;
}

template <typename IndexType, typename MA>
void
switchPackedStorageOrder(StorageOrder orderA, StorageUpLo &upLo,
                         IndexType n, const MA *A, MA *&_A)
{
    if (upLo==Upper) {
        if (orderA==ColMajor) {
            for (IndexType i=0; i<n; ++i) {
                for (IndexType j=i; j<n; ++j) {
                    _A[j+i*(2*n-i-1)/2] = A[i+j*(j+1)/2];
                }
            }
        } else {
            for (IndexType i=0; i<n; ++i) {
                for (IndexType j=i; j<n; ++j) {
                     _A[i+j*(j+1)/2] = A[j+i*(2*n-i-1)/2];
                }
            }
        }
    } else {
        if (orderA==ColMajor) {
            for (IndexType i=0; i<n; ++i) {
                for (IndexType j=0; j<=i; ++j) {
                    _A[j+i*(i+1)/2] = A[i+j*(2*n-j-1)/2];
                }
            }
        } else {
            for (IndexType i=0; i<n; ++i) {
                for (IndexType j=0; j<=i; ++j) {
                    _A[i+j*(2*n-j-1)/2] = A[j+i*(i+1)/2];
                }
            }
        }
    }
}

template <typename IndexType, typename MA>
void
allocateFullStorage(StorageOrder order, IndexType m, IndexType n,
                    MA *&A, IndexType &ldA)
{
    assert(order==RowMajor);

    ldA = n;
    A = new MA[ldA*m];
}

template <typename IndexType, typename MA>
void
allocatePackedStorage(IndexType n, MA *&A)
{
    A = new MA[n*(n+1)/2];
}

template <typename MA>
void
releaseStorage(MA *A)
{
    delete [] A;
}

#endif // CXXBLAS_INTERFACE_AUX_H