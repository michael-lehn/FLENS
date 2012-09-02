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
mv(Transpose trans, const ALPHA &alpha, Permutation P, const DenseVector<VX> &x,
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


} // namespace blas

} // namespace flens
