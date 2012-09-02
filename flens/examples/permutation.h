#include <flens/flens.h>

namespace flens {

struct Permutation
    : public GeneralMatrix<Permutation>
{
    typedef int   ElementType;
    typedef int   IndexType;


    explicit
    Permutation(int n)
        : vector(n)
    {
    }

    DenseVector<Array<int> >  vector;
};


namespace blas {

template <typename ALPHA, typename VX, typename BETA, typename VY>
void
mv(Transpose trans, const ALPHA &alpha,
   const Permutation &P, const DenseVector<VX> &x,
   const BETA &beta, DenseVector<VY> &y)
{
    if (y.length()==0) {
        y.resize(x.length(), x.firstIndex());
    }
    if (trans==NoTrans) {
        int ix0 = x.firstIndex();
        int iy0 = y.firstIndex();
        for (int i=1, ix=ix0, iy=iy0; i<=x.length(); ++i, ++ix, ++iy) {
            y(i) = beta*y(i) + alpha*x(P.vector(i));
        }
    } else {
        int ix0 = x.firstIndex();
        int iy0 = y.firstIndex();
        for (int i=1, ix=ix0, iy=iy0; i<=x.length(); ++i, ++ix, ++iy) {
            y(P.vector(i)) = beta*y(i) + alpha*x(i);
        }
    }
}

template <typename ALPHA, typename MA, typename BETA, typename MB>
void
mm(Transpose transP, Transpose transA, const ALPHA &alpha,
   const Permutation &P, const GeMatrix<MA> &A, const BETA &beta,
   GeMatrix<MB> &B)
{
    typedef typename MB::IndexType   IndexType;

    if (B.numRows()==0 || B.numCols()==0) {
        B.resize(A.numRows(), A.numCols(), A.firstRow(), A.firstCol());
    }

    const Underscore<IndexType> _;

    if (transA==NoTrans) {
        for (IndexType j=1; j<=A.numCols(); ++j) {
            const auto x = A(_,j);
            auto       y = B(_,j);
            mv(transP, alpha, P, x, beta, y);
        }
    } else {
        for (IndexType i=1; i<=A.numRows(); ++i) {
            const auto x = A(i,_);
            auto       y = B(_,i);
            mv(transP, alpha, P, x, beta, y);
        }
    }
}

} // namespace blas

} // namespace flens
