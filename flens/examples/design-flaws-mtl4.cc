#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl; using mtl::iall;

//
// we feed the beast with a const reference
//
void
beast(const dense2D<double> &A)
{
    dense2D<double>   B = sub_matrix(A, 0, 5, 0, 5);
    B[1][1] = 666;
}

int main(int, char**)
{
    using namespace mtl; using mtl::iall;

    const int         m = 5, n = 5;
    dense2D<double>   A(n, n);

    for (int i=0; i<m; ++i) {
        for (int j=0; j<n; ++j) {
            A[i][j] = i+j;
        }
    }

    std::cout << "before: A is\n" << A << "\n";
//
//  calling the beast
//
    beast(A);
    std::cout << "after: A is\n" << A << "\n";

    return 0;
}
