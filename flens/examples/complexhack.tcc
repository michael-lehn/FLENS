#ifndef FLENS_EXAMPLES_COMPLEXHACK_TCC
#define FLENS_EXAMPLES_COMPLEXHACK_TCC 1

#include <flens/examples/complexhack.h>

//
//== RealMatrixClosure =========================================================
//

namespace flens {

template <typename MZ>
RealMatrixClosure<MZ>::RealMatrixClosure(MZ &&Z)
    : _Z(Z)
{
}

template <typename MZ>
template <typename RHS>
void
RealMatrixClosure<MZ>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

template <typename MZ>
RealConstMatrixClosure<MZ>::RealConstMatrixClosure(const MZ &Z)
    : _Z(Z)
{
}

template <typename MZ>
template <typename RHS>
void
RealConstMatrixClosure<MZ>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}


namespace blas {

// RealMatrixClosure -> GeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsGeMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const RealMatrixClosure<MA> &A, MB &&B)
{
    const auto &Z = A._Z;

    ASSERT(trans==NoTrans || trans==Trans);

    if (B.numRows()==0 || B.numCols()==0) {
        if (trans==NoTrans) {
            B.resize(Z.numRows(), Z.numCols(), Z.firstRow(), Z.firstCol());
        } else {
            B.resize(Z.numCols(), Z.numRows(), Z.firstCol(), Z.firstRow());
        }
    }

#   ifndef NDEBUG
    if (trans==NoTrans) {
        ASSERT(Z.numRows()==B.numRows());
        ASSERT(Z.numCols()==B.numCols());
    } else {
        ASSERT(Z.numRows()==B.numCols());
        ASSERT(Z.numCols()==B.numRows());
    }
#   endif

    typedef typename RemoveRef<MB>::Type  MatrixB;
    typedef typename MatrixB::IndexType   IndexType;

    const Underscore<IndexType>  _;

    if (trans==NoTrans) {
        if (B.order()==RowMajor) {
            IndexType  i0 = Z.firstRow();
            IndexType  i1 = Z.lastRow();
            IndexType  I0 = B.firstRow();
            for (IndexType i=i0, I=I0; i<=i1; ++i, ++I) {
                B(I,_) = real(Z(i,_));
            }
        } else {
            IndexType  j0 = Z.firstCol();
            IndexType  j1 = Z.lastCol();
            IndexType  J0 = B.firstCol();
            for (IndexType j=j0, J=J0; j<=j1; ++j, ++J) {
                B(_,J) = real(Z(_,j));
            }
         }
    } else {
        if (B.order()==RowMajor) {
            IndexType  j0 = Z.firstCol();
            IndexType  j1 = Z.lastCol();
            IndexType  I0 = B.firstRow();
            for (IndexType j=j0, I=I0; j<=j1; ++j, ++I) {
                B(I,_) = real(Z(_,j));
            }
       } else {
            IndexType  i0 = Z.firstRow();
            IndexType  i1 = Z.lastRow();
            IndexType  J0 = B.firstCol();
            for (IndexType i=i0, J=J0; i<=i1; ++i, ++J) {
                B(_,J) = real(Z(i,_));
            }
        }
    }
}

// RealConstMatrixClosure -> GeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsGeMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const RealConstMatrixClosure<MA> &A, MB &&B)
{
    std::cerr << "copy2" << std::endl;
}

// GeMatrix -> RealMatrixClosure
template <typename MA, typename MB>
void
copy(Transpose trans, const GeMatrix<MA> &A, RealMatrixClosure<MB> &B)
{
    auto &Z = B._Z;

    ASSERT(trans==NoTrans || trans==Trans);

    if (Z.numRows()==0 || Z.numCols()==0) {
        if (trans==NoTrans) {
            Z.resize(A.numRows(), A.numCols(), A.firstRow(), A.firstCol());
        } else {
            Z.resize(A.numCols(), A.numRows(), A.firstCol(), A.firstRow());
        }
    }

#   ifndef NDEBUG
    if (trans==NoTrans) {
        ASSERT(Z.numRows()==A.numRows());
        ASSERT(Z.numCols()==A.numCols());
    } else {
        ASSERT(Z.numRows()==A.numCols());
        ASSERT(Z.numCols()==A.numRows());
    }
#   endif

    typedef typename MA::IndexType   IndexType;

    const Underscore<IndexType>  _;

    if (trans==NoTrans) {
        if (A.order()==RowMajor) {
            IndexType  i0 = Z.firstRow();
            IndexType  i1 = Z.lastRow();
            IndexType  I0 = A.firstRow();
            for (IndexType i=i0, I=I0; i<=i1; ++i, ++I) {
                real(Z(i,_)) = A(I,_);
            }
        } else {
            IndexType  j0 = Z.firstCol();
            IndexType  j1 = Z.lastCol();
            IndexType  J0 = A.firstCol();
            for (IndexType j=j0, J=J0; j<=j1; ++j, ++J) {
                real(Z(_,j)) = A(_,J);
            }
        }
    } else {
        if (A.order()==RowMajor) {
            IndexType  j0 = Z.firstCol();
            IndexType  j1 = Z.lastCol();
            IndexType  I0 = A.firstRow();
            for (IndexType j=j0, I=I0; j<=j1; ++j, ++I) {
                real(Z(_,j)) = A(I,_);
            }
        } else {
            IndexType  i0 = Z.firstRow();
            IndexType  i1 = Z.lastRow();
            IndexType  J0 = A.firstCol();
            for (IndexType i=i0, J=J0; i<=i1; ++i, ++J) {
                real(Z(i,_)) = A(_,J);
            }
        }
    }
}

} // namespace blas

} // namespace flens


//
//== ImagMatrixClosure =========================================================
//

namespace flens {

template <typename MZ>
ImagMatrixClosure<MZ>::ImagMatrixClosure(MZ &&Z)
    : _Z(Z)
{
}

template <typename MZ>
template <typename RHS>
void
ImagMatrixClosure<MZ>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

template <typename MZ>
ImagConstMatrixClosure<MZ>::ImagConstMatrixClosure(const MZ &Z)
    : _Z(Z)
{
}

template <typename MZ>
template <typename RHS>
void
ImagConstMatrixClosure<MZ>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}


namespace blas {

// ImagMatrixClosure -> GeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsGeMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const ImagMatrixClosure<MA> &A, MB &&B)
{
    const auto &Z = A._Z;

    ASSERT(trans==NoTrans || trans==Trans);

    if (B.numRows()==0 || B.numCols()==0) {
        if (trans==NoTrans) {
            B.resize(Z.numRows(), Z.numCols(), Z.firstRow(), Z.firstCol());
        } else {
            B.resize(Z.numCols(), Z.numRows(), Z.firstCol(), Z.firstRow());
        }
    }

#   ifndef NDEBUG
    if (trans==NoTrans) {
        ASSERT(Z.numRows()==B.numRows());
        ASSERT(Z.numCols()==B.numCols());
    } else {
        ASSERT(Z.numRows()==B.numCols());
        ASSERT(Z.numCols()==B.numRows());
    }
#   endif

    typedef typename RemoveRef<MB>::Type  MatrixB;
    typedef typename MatrixB::IndexType   IndexType;

    const Underscore<IndexType>  _;

    if (trans==NoTrans) {
        if (B.order()==RowMajor) {
            IndexType  i0 = Z.firstRow();
            IndexType  i1 = Z.lastRow();
            IndexType  I0 = B.firstRow();
            for (IndexType i=i0, I=I0; i<=i1; ++i, ++I) {
                B(I,_) = imag(Z(i,_));
            }
        } else {
            IndexType  j0 = Z.firstCol();
            IndexType  j1 = Z.lastCol();
            IndexType  J0 = B.firstCol();
            for (IndexType j=j0, J=J0; j<=j1; ++j, ++J) {
                B(_,J) = imag(Z(_,j));
            }
         }
    } else {
        if (B.order()==RowMajor) {
            IndexType  i0 = Z.firstRow();
            IndexType  i1 = Z.lastRow();
            IndexType  I0 = B.firstCol();
            for (IndexType i=i0, I=I0; i<=i1; ++i, ++I) {
                B(_,I) = imag(Z(i,_));
            }
        } else {
            IndexType  j0 = Z.firstCol();
            IndexType  j1 = Z.lastCol();
            IndexType  J0 = B.firstRow();
            for (IndexType j=j0, J=J0; j<=j1; ++j, ++J) {
                B(J,_) = imag(Z(_,j));
            }
        }
    }
}

// ImagConstMatrixClosure -> GeMatrix
template <typename MA, typename MB>
typename RestrictTo<IsGeMatrix<MB>::value,
         void>::Type
copy(Transpose trans, const ImagConstMatrixClosure<MA> &A, MB &&B)
{
    std::cerr << "copy5" << std::endl;
}

// GeMatrix -> ImagMatrixClosure
template <typename MA, typename MB>
void
copy(Transpose trans, const GeMatrix<MA> &A, ImagMatrixClosure<MB> &B)
{
    auto &Z = B._Z;

    ASSERT(trans==NoTrans || trans==Trans);

    if (Z.numRows()==0 || Z.numCols()==0) {
        if (trans==NoTrans) {
            Z.resize(A.numRows(), A.numCols(), A.firstRow(), A.firstCol());
        } else {
            Z.resize(A.numCols(), A.numRows(), A.firstCol(), A.firstRow());
        }
    }

#   ifndef NDEBUG
    if (trans==NoTrans) {
        ASSERT(Z.numRows()==A.numRows());
        ASSERT(Z.numCols()==A.numCols());
    } else {
        ASSERT(Z.numRows()==A.numCols());
        ASSERT(Z.numCols()==A.numRows());
    }
#   endif

    typedef typename MA::IndexType   IndexType;

    const Underscore<IndexType>  _;

    if (trans==NoTrans) {
        if (A.order()==RowMajor) {
            IndexType  i0 = Z.firstRow();
            IndexType  i1 = Z.lastRow();
            IndexType  I0 = A.firstRow();
            for (IndexType i=i0, I=I0; i<=i1; ++i, ++I) {
                imag(Z(i,_)) = A(I,_);
            }
        } else {
            IndexType  j0 = Z.firstCol();
            IndexType  j1 = Z.lastCol();
            IndexType  J0 = A.firstCol();
            for (IndexType j=j0, J=J0; j<=j1; ++j, ++J) {
                imag(Z(_,j)) = A(_,J);
            }
        }
    } else {
        if (A.order()==RowMajor) {
            IndexType  j0 = Z.firstCol();
            IndexType  j1 = Z.lastCol();
            IndexType  I0 = A.firstRow();
            for (IndexType j=j0, I=I0; j<=j1; ++j, ++I) {
                imag(Z(_,j)) = A(I,_);
            }
        } else {
            IndexType  i0 = Z.firstRow();
            IndexType  i1 = Z.lastRow();
            IndexType  J0 = A.firstCol();
            for (IndexType i=i0, J=J0; i<=i1; ++i, ++J) {
                imag(Z(i,_)) = A(_,J);
            }
        }
    }
}

} // namespace blas

} // namespace flens


//
//== real ======================================================================
//

namespace flens {

template <typename VZ>
typename RealVector<DenseVector<VZ> >::ConstView
real(const DenseVector<VZ> &z)
{
    typedef typename RealVector<DenseVector<VZ> >::ConstViewEngine  Engine;
    typedef typename VZ::ElementType                                T;
    typedef typename Engine::ElementType                            PT;

    static_assert(sizeof(T)%sizeof(PT)==0, "Are you kidding me?");

    return Engine(z.length(),
                  reinterpret_cast<const PT *>(z.data()),
                  z.stride()*sizeof(T)/sizeof(PT),
                  z.firstIndex());
}

template <typename VZ>
typename RestrictTo<IsDenseVector<VZ>::value,
         typename RealVector<typename RemoveRef<VZ>::Type>::View>::Type
real(VZ &&z)
{
    typedef typename RemoveRef<VZ>::Type                VectorZ;
    typedef typename RealVector<VectorZ>::ViewEngine    Engine;
    typedef typename VectorZ::ElementType               T;
    typedef typename Engine::ElementType                PT;

    static_assert(sizeof(T)%sizeof(PT)==0, "Are you kidding me?");

    return Engine(z.length(),
                  reinterpret_cast<PT *>(z.data()),
                  z.stride()*sizeof(T)/sizeof(PT),
                  z.firstIndex());
}

template <typename MZ>
RealConstMatrixClosure<MZ>
real(const GeMatrix<MZ> &Z)
{
    return Z;
}

template <typename MZ>
typename RestrictTo<IsGeMatrix<MZ>::value,
         RealMatrixClosure<MZ> >::Type
real(MZ &&Z)
{
    return Z;
}

} // namespace flens


//
//== imag ======================================================================
//

namespace flens {

template <typename VZ>
typename ImagVector<DenseVector<VZ> >::ConstView
imag(const DenseVector<VZ> &z)
{
    typedef typename ImagVector<DenseVector<VZ> >::ConstViewEngine  Engine;
    typedef typename VZ::ElementType                                T;
    typedef typename Engine::ElementType                            PT;

    static_assert(sizeof(T)==2*sizeof(PT), "Are you kidding me?");

    return Engine(z.length(),
                  reinterpret_cast<const PT *>(z.data()) + 1,
                  2*z.stride(),
                  z.firstIndex());
}

template <typename VZ>
typename RestrictTo<IsDenseVector<VZ>::value,
         typename ImagVector<typename RemoveRef<VZ>::Type>::View>::Type
imag(VZ &&z)
{
    typedef typename RemoveRef<VZ>::Type                VectorZ;
    typedef typename ImagVector<VectorZ>::ViewEngine    Engine;
    typedef typename VectorZ::ElementType               T;
    typedef typename Engine::ElementType                PT;

    static_assert(sizeof(T)%sizeof(PT)==0, "Are you kidding me?");

    return Engine(z.length(),
                  reinterpret_cast<PT *>(z.data()) + 1,
                  2*z.stride(),
                  z.firstIndex());
}

template <typename MZ>
ImagConstMatrixClosure<MZ>
imag(const GeMatrix<MZ> &Z)
{
    return Z;
}

template <typename MZ>
typename RestrictTo<IsGeMatrix<MZ>::value,
         ImagMatrixClosure<MZ> >::Type
imag(MZ &&Z)
{
    return Z;
}


} // namespace flens

#endif // FLENS_EXAMPLES_COMPLEXHACK_TCC
