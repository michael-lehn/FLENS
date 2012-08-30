#ifndef FLENS_EXAMPLES_COMPLEXHACK_H
#define FLENS_EXAMPLES_COMPLEXHACK_H 1

#include <flens/flens.h>

//
//== RealArray =================================================================
//

namespace flens {

template <typename A>
struct RealArray
{
};

template <typename T, typename I, typename A>
struct RealArray<Array<T, I, A> >
{
    typedef typename ComplexTrait<T>::PrimitiveType   PT;
    typedef ConstArrayView<PT, I, A>                  ConstView;
    typedef ArrayView<PT, I, A>                       View;
};

template <typename T, typename I, typename A>
struct RealArray<ArrayView<T, I, A> >
{
    typedef typename ComplexTrait<T>::PrimitiveType   PT;
    typedef ConstArrayView<PT, I, A>                  ConstView;
    typedef ArrayView<PT, I, A>                       View;
};

template <typename T, typename I, typename A>
struct RealArray<ConstArrayView<T, I, A> >
{
    typedef typename ComplexTrait<T>::PrimitiveType   PT;
    typedef ConstArrayView<PT, I, A>                  ConstView;
    typedef ArrayView<PT, I, A>                       View;
};

} // namespace flens


//
//== ImagArray =================================================================
//

namespace flens {

template <typename A>
struct ImagArray
{
};

template <typename PT, typename I, typename A>
struct ImagArray<Array<std::complex<PT>, I, A> >
{
    typedef ConstArrayView<PT, I, A>                  ConstView;
    typedef ArrayView<PT, I, A>                       View;
};

template <typename PT, typename I, typename A>
struct ImagArray<ArrayView<std::complex<PT>, I, A> >
{
    typedef ConstArrayView<PT, I, A>                  ConstView;
    typedef ArrayView<PT, I, A>                       View;
};

template <typename PT, typename I, typename A>
struct ImagArray<ConstArrayView<std::complex<PT>, I, A> >
{
    typedef ConstArrayView<PT, I, A>                  ConstView;
    typedef ArrayView<PT, I, A>                       View;
};

} // namespace flens


//
//== RealVector ================================================================
//

namespace flens {

template <typename VZ>
struct RealVector
{
};

template <typename VZ>
struct RealVector<DenseVector<VZ> >
{
    typedef typename RealArray<VZ>::ConstView   ConstViewEngine;
    typedef typename RealArray<VZ>::View        ViewEngine;
    typedef DenseVector<ConstViewEngine>        ConstView;
    typedef DenseVector<ViewEngine>             View;
};

} // namespace flens


//
//== ImagVector ================================================================
//

namespace flens {

template <typename VZ>
struct ImagVector
{
};

template <typename VZ>
struct ImagVector<DenseVector<VZ> >
{
    typedef typename ImagArray<VZ>::ConstView   ConstViewEngine;
    typedef typename ImagArray<VZ>::View        ViewEngine;
    typedef DenseVector<ConstViewEngine>        ConstView;
    typedef DenseVector<ViewEngine>             View;
};

} // namespace flens


//
//== RealMatrixClosure =========================================================
//

namespace flens {

template <typename MZ>
struct RealMatrixClosure
    : public GeneralMatrix<RealMatrixClosure<MZ> >
{
    typedef typename RemoveRef<MZ>::Type    MatrixZ;
    typedef typename MatrixZ::ElementType   ElementType;
    typedef typename MatrixZ::IndexType     IndexType;

    RealMatrixClosure(MZ &&Z);

    template <typename RHS>
        void
        operator=(const Matrix<RHS> &rhs);

    MZ &&_Z;
};

template <typename MZ>
struct RealConstMatrixClosure
    : public GeneralMatrix<RealConstMatrixClosure<MZ> >
{
    typedef typename MZ::ElementType   ElementType;
    typedef typename MZ::IndexType     IndexType;

    RealConstMatrixClosure(const MZ &Z);

    template <typename RHS>
        void
        operator=(const Matrix<RHS> &rhs);

    const MZ &_Z;
};

namespace blas {

// RealMatrixClosure -> GeMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsGeMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const RealMatrixClosure<MA> &A, MB &&B);

// RealConstMatrixClosure -> GeMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsGeMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const RealConstMatrixClosure<MA> &A, MB &&B);

// GeMatrix -> RealMatrixClosure
template <typename MA, typename MB>
    void
    copy(Transpose trans, const GeMatrix<MA> &A, RealMatrixClosure<MB> &B);

} // namespace blas

} // namespace flens


//
//== ImagMatrixClosure =========================================================
//

namespace flens {

template <typename MZ>
struct ImagMatrixClosure
    : public GeneralMatrix<ImagMatrixClosure<MZ> >
{
    typedef typename RemoveRef<MZ>::Type    MatrixZ;
    typedef typename MatrixZ::ElementType   ElementType;
    typedef typename MatrixZ::IndexType     IndexType;

    ImagMatrixClosure(MZ &&Z);

    template <typename RHS>
        void
        operator=(const Matrix<RHS> &rhs);

    MZ &&_Z;
};

template <typename MZ>
struct ImagConstMatrixClosure
    : public GeneralMatrix<ImagConstMatrixClosure<MZ> >
{
    typedef typename MZ::ElementType   ElementType;
    typedef typename MZ::IndexType     IndexType;

    ImagConstMatrixClosure(const MZ &Z);

    template <typename RHS>
        void
        operator=(const Matrix<RHS> &rhs);

    const MZ &_Z;
};


namespace blas {

// ImagMatrixClosure -> GeMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsGeMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const ImagMatrixClosure<MA> &A, MB &&B);

// ImagConstMatrixClosure -> GeMatrix
template <typename MA, typename MB>
    typename RestrictTo<IsGeMatrix<MB>::value,
             void>::Type
    copy(Transpose trans, const ImagConstMatrixClosure<MA> &A, MB &&B);

// GeMatrix -> ImagMatrixClosure
template <typename MA, typename MB>
    void
    copy(Transpose trans, const GeMatrix<MA> &A, ImagMatrixClosure<MB> &B);

} // namespace blas

} // namespace flens


//
//== real ======================================================================
//

namespace flens {

template <typename VZ>
    typename RealVector<DenseVector<VZ> >::ConstView
    real(const DenseVector<VZ> &z);

template <typename VZ>
    typename RestrictTo<IsDenseVector<VZ>::value,
             typename RealVector<typename RemoveRef<VZ>::Type>::View>::Type
    real(VZ &&z);

template <typename MZ>
    RealConstMatrixClosure<MZ>
    real(const GeMatrix<MZ> &Z);

template <typename MZ>
    typename RestrictTo<IsGeMatrix<MZ>::value,
             RealMatrixClosure<MZ> >::Type
    real(MZ &&Z);

} // namespace flens


//
//== imag ======================================================================
//

namespace flens {

template <typename VZ>
    typename ImagVector<DenseVector<VZ> >::ConstView
    imag(const DenseVector<VZ> &z);

template <typename VZ>
    typename RestrictTo<IsDenseVector<VZ>::value,
             typename ImagVector<typename RemoveRef<VZ>::Type>::View>::Type
    imag(VZ &&z);

template <typename MZ>
    ImagConstMatrixClosure<MZ>
    imag(const GeMatrix<MZ> &Z);

template <typename MZ>
    typename RestrictTo<IsGeMatrix<MZ>::value,
             ImagMatrixClosure<MZ> >::Type
    imag(MZ &&Z);

} // namespace flens

#endif // FLENS_EXAMPLES_COMPLEXHACK_H
