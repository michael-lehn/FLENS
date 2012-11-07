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
template <typename X>
struct Eval
{
    template <typename IndexType>
    static typename X::ElementType
    element(const X &x, IndexType i)
    {
        return x(i);
    }
};

template <typename VX, typename VY>
struct Eval<VectorClosure<OpAdd, VX, VY> >
{
    template <typename IndexType>
    static typename VectorClosure<OpAdd, VX, VY>::ElementType
    element(const VectorClosure<OpAdd, VX, VY> &x, IndexType i)
    {
        return Eval<VX>::element(x.left(), i)
             + Eval<VY>::element(x.right(), i);
    }
};

template <typename SV, typename VX>
struct Eval<VectorClosure<OpMult, SV, VX> >
{
    template <typename IndexType>
    static typename VectorClosure<OpMult, SV, VX>::ElementType
    element(const VectorClosure<OpMult, SV, VX> &x, IndexType i)
    {
        return x.left().value()*Eval<VX>::element(x.right(), i);
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
        y(i) = Eval<VC>::element(x, i);
    }
}

//
//  Entry point for the plus-assignment of vector sums detected by
//  `IsDenseVectorSum`
//
template <typename ALPHA, typename VC, typename VY>
typename RestrictTo<IsDenseVectorSum<VC>::value
                 && IsDenseVector<VY>::value,
         void>::Type
axpy(const ALPHA &alpha, const VC &x, VY &&y)
{
    typedef typename RemoveRef<VY>::Type  VectorY;
    typedef typename VectorY::IndexType   IndexType;

    for (IndexType i=y.firstIndex(); i<=y.lastIndex(); ++i) {
        y(i) += alpha*Eval<VC>::element(x, i);
    }
}

} } // namespace blas, flens

#endif // FLENS_EXAMPLES_SIMPLE_ET_H
