#ifndef FLENS_EXAMPLES_SIMPLE_ET_H
#define FLENS_EXAMPLES_SIMPLE_ET_H 1

#include <flens/flens.h>

//------------------------------------------------------------------------------

namespace flens { namespace blas {

//
//  This trait detects whether a vector closures encapsulates vector sums of
//  the form alpha1*x1 + alpha2*x2 + alpha3*x3 ...
//
template <typename V>
struct IsDenseVectorSum
{
    static const bool value = IsDenseVector<V>::value;
};

template <typename L, typename R>
struct IsDenseVectorSum<VectorClosure<OpAdd, L, R> >
{
    static const bool value = IsDenseVectorSum<L>::value
                           && IsDenseVectorSum<R>::value;
};

template <typename SV, typename VX>
struct IsDenseVectorSum<VectorClosure<OpMult, SV, VX> >
{
    static const bool value = IsScalarValue<SV>::value
                           && IsDenseVectorSum<VX>::value;
};

//
//  This trait defines that we want to use a non-default evaluation for
//  vector sums detected by the `IsDenseVectorSum` trait.
//
template <typename Op, typename L, typename R>
struct DefaultEval<VectorClosure<Op, L, R> >
{
    typedef VectorClosure<OpAdd, L, R>   Closure;

    static const bool value = !IsDenseVectorSum<Closure>::value;
};

//
//  Trait for the component-wise evaluation of closures
//
template <typename VX, typename VY>
struct Eval
{
    template <typename IndexType>
    static void
    copy(IndexType i, const VX &x, VY &y)
    {
        y(i) = x(i);
    }

    template <typename IndexType, typename ALPHA>
    static void
    axpy(IndexType i, const ALPHA &alpha, const VX &x, VY &y)
    {
        y(i) += alpha*x(i);
    }
};

template <typename VL, typename VR, typename VY>
struct Eval<VectorClosure<OpAdd, VL, VR>, VY>
{
    template <typename IndexType>
    static void
    copy(IndexType i, const VectorClosure<OpAdd, VL, VR> &x, VY &y)
    {
        typedef typename VY::ElementType  ElementType;
        const ElementType One(1);

        Eval<VL, VY>::copy(i, x.left(), y);
        Eval<VR, VY>::axpy(i, One, x.right(), y);
    }

    template <typename IndexType, typename ALPHA>
    static void
    axpy(IndexType i, const ALPHA &alpha, const VectorClosure<OpAdd, VL, VR> &x,
         VY &y)
    {
        typedef typename VY::ElementType  ElementType;
        const ElementType One(1);

        Eval<VL, VY>::axpy(i, alpha, x.left(), y);
        Eval<VR, VY>::axpy(i, alpha, x.right(), y);
    }
};

template <typename SV, typename VX, typename VY>
struct Eval<VectorClosure<OpMult, SV, VX>, VY>
{
    template <typename IndexType>
    static void
    copy(IndexType i, const VectorClosure<OpMult, SV, VX> &x, VY &y)
    {
        Eval<VX, VY>::copy(i, x.right(), y);
        y(i) *= x.left().value();
    }

    template <typename IndexType, typename ALPHA>
    static void
    axpy(IndexType i, const ALPHA &alpha,
         const VectorClosure<OpMult, SV, VX> &x, VY &y)
    {
        Eval<VX, VY>::axpy(i, alpha*x.left().value(), x.right(), y);
    }
};

//
//  Entry point for the assignment of vector sums detected by `IsDenseVectorSum`
//
template <typename VC, typename VY>
typename RestrictTo<IsDenseVectorSum<VC>::value
                 && IsDenseVector<VY>::value,
         void>::Type
copy(const VC &x, VY &&y)
{
    typedef typename RemoveRef<VY>::Type  VectorY;
    typedef typename VectorY::IndexType   IndexType;

    for (IndexType i=y.firstIndex(); i<=y.lastIndex(); ++i) {
        Eval<VC, VectorY>::copy(i, x, y);
    }
}

template <typename ALPHA, typename VC, typename VY>
typename RestrictTo<IsDenseVectorSum<VC>::value
                 && IsDenseVector<VY>::value,
         void>::Type
axpy(const ALPHA &alpha, const VC &x, VY &&y)
{
    typedef typename RemoveRef<VY>::Type  VectorY;
    typedef typename VectorY::IndexType   IndexType;

    for (IndexType i=y.firstIndex(); i<=y.lastIndex(); ++i) {
        Eval<VC, VectorY>::axpy(i, alpha, x, y);
    }
}

} } // namespace blas, flens

#endif // FLENS_EXAMPLES_SIMPLE_ET_H
