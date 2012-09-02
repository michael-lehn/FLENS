#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    const int n = 10;
    typedef TinyArray<double, n>                Engine;
    typedef TinyArrayView<double, n/2, 2>       EngineView;
    typedef TinyConstArrayView<double, n/2, 2>  EngineConstView;
    typedef TinyVector<Engine>                  TVector;
    typedef TinyVector<EngineView>              TVectorView;
    typedef TinyVector<EngineConstView>         TVectorConstView;

    TVector x, y;

    for (int i=1; i<=n; ++i) {
        x(i) = i;
        y(i) = 0;
    }

    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    y += double(2)*x;
    cout << "y = " << y << endl;

    const TVectorView  x2(EngineView(x.data()));
    TVectorView        y2(EngineView(y.data()));

    cout << "x2 = " << x2 << endl;
    cout << "y2 = " << y2 << endl;
    y2 += double(2)*x2;
    cout << "y = " << y << endl;

    TVectorConstView  cx2(EngineConstView(x.engine().data()));
    cout << "cx2 = " << cx2 << endl;

    y2 += cx2;
    cout << "y = " << y << endl;

    y *= 2;
    cout << "y = " << y << endl;

    y /= 2;
    cout << "y = " << y << endl;

    y = 2*x + 4*x + 7*x;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
}
