#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    const int n = 20;

    DenseVector<Array<double> >  x(n);

    Underscore<int> _;

    for (int i=1; i<=x.length(); ++i) {
        x(i) = i;
    }

    cout << "x = " << x << endl;

    GeMatrix<FullStorageView<double, ColMajor> >  A1(4, 5, x);
    GeMatrix<FullStorageView<double, RowMajor> >  A2(4, 5, x);

    cout << "A1 = " << A1 << endl;
    cout << "A2 = " << A2 << endl;


    GeMatrix<FullStorageView<double, ColMajor> >  B1(2, 5, x(_(1,2,n)));
    GeMatrix<FullStorageView<double, RowMajor> >  B2(2, 5, x(_(1,2,n)));

    cout << "B1 = " << B1 << endl;
    cout << "B2 = " << B2 << endl;

    GeMatrix<FullStorageView<double, ColMajor> >  C1(2, 3, x(_(1,2,n)), 2, 2);
    GeMatrix<FullStorageView<double, RowMajor> >  C2(2, 3, x(_(1,2,n)), 2, 2);

    cout << "C1 = " << C1 << endl;
    cout << "C2 = " << C2 << endl;



}
