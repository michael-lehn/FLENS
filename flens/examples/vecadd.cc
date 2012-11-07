#include <iostream>
#include <flens/flens.h>

#ifdef VECSUM

#define NB 8

//------------------------------------------------------------------------------

namespace flens {

template <typename ALPHA, typename TX, typename TY>
void
copy_NB(ALPHA alpha, const TX *x, TY *y)
{
    /*
    for (int i=0; i<NB; ++i) {
        y[i] = alpha*x[i];
    }
    */
    y[ 0] = x[ 0];
    y[ 1] = x[ 1];
    y[ 2] = x[ 2];
    y[ 3] = x[ 3];
    y[ 4] = x[ 4];
    y[ 5] = x[ 5];
    y[ 6] = x[ 6];
    y[ 7] = x[ 7];

    y[ 0] *= alpha;
    y[ 1] *= alpha;
    y[ 2] *= alpha;
    y[ 3] *= alpha;
    y[ 4] *= alpha;
    y[ 5] *= alpha;
    y[ 6] *= alpha;
    y[ 7] *= alpha;
}

template <typename ALPHA, typename TX, typename TY>
void
add_NB(ALPHA alpha, const TX *x, TY *y)
{
    /*
    for (int i=0; i<NB; ++i) {
        y[i] += alpha*x[i];
    }
    */

#   ifndef HAVE_SSE
    y[ 0] += alpha*x[ 0];
    y[ 1] += alpha*x[ 1];
    y[ 2] += alpha*x[ 2];
    y[ 3] += alpha*x[ 3];
    y[ 4] += alpha*x[ 4];
    y[ 5] += alpha*x[ 5];
    y[ 6] += alpha*x[ 6];
    y[ 7] += alpha*x[ 7];
#   else
#       ifdef HAVE_AVX

        __m256d x1_, x2_, y1_, y2_;
        __m256d alpha_ = _mm256_set1_pd(alpha);

        x1_ = _mm256_loadu_pd(x);
        x2_ = _mm256_loadu_pd(x+4);

        x1_ = _mm256_mul_pd(alpha_, x1_);
        x2_ = _mm256_mul_pd(alpha_, x2_);

        y1_ = _mm256_loadu_pd(y);
        y2_ = _mm256_loadu_pd(y+4);

        y1_ = _mm256_add_pd(y1_, x1_);
        y2_ = _mm256_add_pd(y2_, x2_);

        _mm256_storeu_pd(y  , y1_);
        _mm256_storeu_pd(y+4, y2_);   

#       else

        __m128d x1_, x2_, x3_, x4_ , y1_, y2_, y3_, y4_;
        __m128d alpha_ = _mm_set1_pd(alpha);

        x1_ = _mm_loadu_pd(x);
        x2_ = _mm_loadu_pd(x+2);
        x3_ = _mm_loadu_pd(x+4);
        x4_ = _mm_loadu_pd(x+6);

        x1_ = _mm_mul_pd(alpha_, x1_);
        x2_ = _mm_mul_pd(alpha_, x2_);
        x3_ = _mm_mul_pd(alpha_, x3_);
        x4_ = _mm_mul_pd(alpha_, x4_);

        y1_ = _mm_loadu_pd(y);
        y2_ = _mm_loadu_pd(y+2);
        y3_ = _mm_loadu_pd(y+4);
        y4_ = _mm_loadu_pd(y+6);

        y1_ = _mm_add_pd(y1_, x1_);
        y2_ = _mm_add_pd(y2_, x2_);
        y3_ = _mm_add_pd(y3_, x3_);
        y4_ = _mm_add_pd(y4_, x4_);

        _mm_storeu_pd(y  , y1_);
        _mm_storeu_pd(y+2, y2_);
        _mm_storeu_pd(y+4, y3_);
        _mm_storeu_pd(y+6, y4_);  

#       endif // HAVE_AVX
#   endif // HAVE_SSE

}

} // namespace flens

//------------------------------------------------------------------------------

namespace flens {

template <typename IndexType, typename ALPHA, typename VX, typename VY>
void
copy(IndexType i0, IndexType nb, ALPHA alpha, const DenseVector<VX> &x,
     DenseVector<VY> &y)
{
    typedef typename DenseVector<VX>::ElementType  TX;
    typedef typename DenseVector<VY>::ElementType  TY;

    const IndexType  incX = x.stride();
    const TX         *_x = x.data() + i0*incX;

    const IndexType  incY = y.stride();
    TY               *_y  = y.data() + i0*incY;

    if (incX==1 && incY==1) {
        if (nb==NB) {
            copy_NB(alpha, _x, _y);
        } else {
            for (IndexType i=0; i<nb; ++i) {
                _y[i] = alpha * _x[i];
            }
        }
    } else {
        for (IndexType i=0; i<nb; ++i, _x+=incX, _y+=incY) {
            *_y = alpha * (*_x);
        }
    }
}

template <typename IndexType, typename ALPHA, typename VX, typename VY>
void
add(IndexType i0, IndexType nb, ALPHA alpha, const DenseVector<VX> &x,
    DenseVector<VY> &y)
{
    typedef typename DenseVector<VX>::ElementType  TX;
    typedef typename DenseVector<VY>::ElementType  TY;

    const IndexType  incX = x.stride();
    const TX         *_x = x.data() + i0*incX;

    const IndexType  incY = y.stride();
    TY               *_y  = y.data() + i0*incY;

    if (incX==1 && incY==1) {
        if (nb==NB) {
            add_NB(alpha, _x, _y);
        } else {
            for (IndexType i=0; i<nb; ++i) {
                _y[i] += alpha * _x[i];
            }
        }
    } else {
        for (IndexType i=0; i<nb; ++i, _x+=incX, _y+=incY) {
            *_y += alpha * (*_x);
        }
    }
}

} // namespace flens

//------------------------------------------------------------------------------

namespace flens {

template <typename IndexType, typename ALPHA, typename VX1, typename VX2,
          typename VY>
void
copy(IndexType i0, IndexType nb, ALPHA alpha,
     const VectorClosure<OpAdd, VX1, VX2> &x,
     DenseVector<VY> &y)
{
    copy(i0, nb, alpha, x.left(), y);
    add(i0, nb, alpha, x.right(), y);
}

template <typename IndexType, typename ALPHA, typename SV, typename VX,
          typename VY>
typename RestrictTo<IsScalarValue<SV>::value,
                    void>::Type
copy(IndexType i0, IndexType nb, ALPHA alpha,
     const VectorClosure<OpMult, SV, VX> &x,
     DenseVector<VY> &y)
{
    copy(i0, nb, alpha*x.left().value(), x.right(), y);
}

template <typename IndexType, typename ALPHA, typename VX1, typename VX2,
          typename VY>
void
add(IndexType i0, IndexType nb, ALPHA alpha,
    const VectorClosure<OpAdd, VX1, VX2> &x,
    DenseVector<VY> &y)
{
    add(i0, nb, alpha, x.left(), y);
    add(i0, nb, alpha, x.right(), y);
}

template <typename IndexType, typename ALPHA, typename VX1, typename VX2,
          typename VY>
void
add(IndexType i0, IndexType nb, ALPHA alpha,
    const VectorClosure<OpMult, VX1, VX2> &x,
    DenseVector<VY> &y)
{
    add(i0, nb, alpha*x.left().value(), x.right(), y);
}

} // namespace flens

//------------------------------------------------------------------------------

namespace flens { namespace blas {

//
//  IsDenseVectorSum: Trait for vector sums of form alpha1*x1 + alpha2*x2 + ...
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
//  Use non-default evaluation for vector sums
//
template <typename Op, typename L, typename R>
struct DefaultEval<VectorClosure<Op, L, R> >
{
    typedef VectorClosure<OpAdd, L, R>   Closure;

    static const bool value = !IsDenseVectorSum<Closure>::value;
};

//
//  Entry point for the evaluation of vector sums.
//
template <typename VL, typename VR, typename VY>
typename RestrictTo<IsDenseVectorSum<VectorClosure<OpAdd, VL, VR> >::value,
         void>::Type
copy(const VectorClosure<OpAdd, VL, VR> &x, DenseVector<VY> &y)
{
    typedef typename DenseVector<VY>::IndexType  IndexType;
    typedef typename DenseVector<VY>::ElementType  ElementType;

    const IndexType nb = NB;
    const IndexType n  = y.length();
    const IndexType n1 = (n/nb)*nb;

    const ElementType  One(1);

    for (IndexType i=0; i<n1; i+=nb) {
        copy(i, nb, One, x.left(), y);
        add(i, nb, One, x.right(), y);
    }
    if (n1<n) {
        copy(n1, n-n1, One, x.left(), y);
        add(n1, n-n1, One, x.right(), y);
    }
}

} } // namespace blas, flens

#endif

//------------------------------------------------------------------------------

namespace flens {

template <typename IndexType,
          typename ALPHA1, typename VX1,
          typename ALPHA2, typename VX2,
          typename ALPHA3, typename VX3,
          typename VY>
void
sum(IndexType n,
    const ALPHA1 alpha1, const VX1 *x1, IndexType incX1,
    const ALPHA2 alpha2, const VX2 *x2, IndexType incX2,
    const ALPHA3 alpha3, const VX3 *x3, IndexType incX3,
    VY *y, IndexType incY)
{
    for (IndexType i=0; i<n; ++i, x1+=incX1, x2+=incX2, x3+=incX3, y+=incY) {
        *y = alpha1 * (*x1) + alpha2 * (*x2) + alpha3 * (*x3);
    }
}


template <typename ALPHA1, typename VX1,
          typename ALPHA2, typename VX2,
          typename ALPHA3, typename VX3,
          typename VY>
void
sum(const ALPHA1 alpha1, const DenseVector<VX1> &x1,
    const ALPHA2 alpha2, const DenseVector<VX2> &x2,
    const ALPHA3 alpha3, const DenseVector<VX3> &x3,
    DenseVector<VY> &y)
{
    sum(y.length(),
        alpha1, x1.data(), x1.stride(),
        alpha2, x2.data(), x2.stride(),
        alpha3, x3.data(), x3.stride(),
        y.data(), y.stride());
}

} // namespace flens

//------------------------------------------------------------------------------

#include <iostream>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    const int n = 150;
    DenseVector<Array<double> >   x1(n), x2(n), x3(n), y(n);

    x1 = 1;
    x2 = 2;
    x3 = 3;

    for (int count=1; count<=500000; ++count) {

        sum(2.0, x1, 4.0, x2, 3.0, x3, y);
        sum(2.0, x1, 4.0, x2, 3.0, x3, y);
        sum(2.0, x1, 4.0, x2, 3.0, x3, y);
        sum(2.0, x1, 4.0, x2, 3.0, x3, y);

        /*
        y = 2.0*x1 + 4.0*x2 + 3.0*x3;
        y = 2.0*x1 + 4.0*x2 + 3.0*x3;
        y = 2.0*x1 + 4.0*x2 + 3.0*x3;
        y = 2.0*x1 + 4.0*x2 + 3.0*x3;
        */

    }

    return y(2);
}
