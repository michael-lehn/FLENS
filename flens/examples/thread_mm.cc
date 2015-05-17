#include <flens/flens.cxx>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>

using namespace flens;
using namespace std;

void
mm(const GeMatrix<FullStorage<double> > *A,
   const GeMatrix<FullStorageView<double, ColMajor> > *B,
   GeMatrix<FullStorageView<double, ColMajor> > *C)
{
    blas::mm(NoTrans, NoTrans, 1.0, *A, *B, 0.0, *C);
}

int
main()
{
    const int m = 4000;
    const int k = 4000;
    const int n = 4000;

    GeMatrix<FullStorage<double> > A(m, k), B(k,n), C(m,n);

    const Underscore<int>  _;

    fillRandom(A);
    fillRandom(B);

    auto C_1 = C(_, _(    1, n/2));
    auto C_2 = C(_, _(n/2+1,   n));

    const auto B_1 = B(_, _(    1, n/2));
    const auto B_2 = B(_, _(n/2+1,   n));

#   ifndef USE_THREADS
    C = A*B;
#   else

    std::vector<std::function<void()> >   task;

    task.push_back([=, &A]() mutable { C_1 = A*B_1; });
    task.push_back([=, &A]() mutable { C_2 = A*B_2; });

    auto p0 = std::thread(task[0]);
    auto p1 = std::thread(task[1]);

    p0.join();
    p1.join();
#   endif

    cout << "asum(A*B) = " << blas::asum(C) << endl;
}
