#include <flens/flens.h>

namespace flens {


///
/// This is a minimalistic user-defined matrix type:
///  - It's a *symmetric type* as we derive it from *SymmetricMatrix*
///  - A FLENS matrix type needs to define at least the two typedefs
///    *IndexType* and *ElementType*
///  - Otherwise this matrix type is pretty useless so far.
///
struct MySyMatrix
    : public SymmetricMatrix<MySyMatrix>
{
    typedef double ElementType;
    typedef int    IndexType;
};

} // namespace flens


namespace flens { namespace blas {

///
/// Now we give it some meaning.  We define a matrix-vector product for it
/// where the vector is of type *DenseVector*.  The Operation we define has
/// the form $y \leftarrow \beta y + \alpha A x$.
///
/// *We define this product as if the matrix $A$ has the special tridiagonal
/// form:*
///
///  *--[LATEX]---------------------------------------------*
///  |                                                      |
///  |   \begin{pmatrix}                                    |
///  |    2  & -1     & 0      &        & 0      &  0 \\    |
///  |   -1  & \ddots & \ddots & \ddots &        &  0 \\    |
///  |    0  & \ddots & \ddots & \ddots & \ddots &    \\    |
///  |       & \ddots & \ddots & \ddots & \ddots &  0 \\    |
///  |    0  &        & \ddots & \ddots & \ddots & -1 \\    |
///  |    0  & 0      &        &      0 & -1     &  2       |
///  |   \end{pmatrix}                                      |
///  |                                                      |
///  *------------------------------------------------------*
///
/// The matrix-vector product gets computed by this function.
///
template <typename ALPHA, typename VX, typename BETA, typename VY>
void
mv(const ALPHA &alpha, const MySyMatrix &, const DenseVector<VX> &x,
   const BETA &beta, DenseVector<VY> &y)
///-
{
    typedef typename DenseVector<VX>::IndexType  IndexType;
    const IndexType n = x.length();

    // we only resize empty left hand sides
    if (y.length()==0) {
        y.resize(n);
    }
    ASSERT(y.length()==n);


    for (IndexType i=1; i<=n; ++i) {
        y(i) = beta*y(i) + alpha*2*x(i);
        if (i>1) {
            y(i) -= alpha*x(i-1);
        }
        if (i<n) {
            y(i) -= alpha*x(i+1);
        }
    }
}

} } // namespace blas, flens
