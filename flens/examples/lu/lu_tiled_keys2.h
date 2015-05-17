#ifndef LU_LU_TILED_KEYS2_H
#define LU_LU_TILED_KEYS2_H 1

#include <flens/examples/lu/format.h>

namespace flens {

// range: 0 <= key < mn = min(mn)
auto keyLU      = [=](int i)
                  {
                      return i-1;
                  };

// range: min(m,n) <= key < min(mn,n) + m*n
auto keyPtA     = [=](int mn, int m, int i, int j)
                  {
                      return m*(j-1)+(i-1)+mn;
                  };

// range: min(mn,n) + m*n <= key < 2*min(mn,n) +m*n
auto keyUpdateP = [=](int mn, int m, int n, int i)
                  {
                      return (i-1)+mn+m*n;
                  };

// range: 2*min(mn,n) +m*n <= key < 2*min(mn,n) +2*m*n
auto keyApplyU  = [=](int mn, int m, int n, int i, int j)
                  {
                      return m*(j-1)+(i-1) + 2*mn + m*n;
                  };

// range: 2*min(mn,n) +2m*n <= key < 2*min(mn,n) +3*m*n
auto keyApplyL  = [=](int mn, int m, int n, int i, int j)
                  {
                      return m*(j-1)+(i-1) + 2*mn + 2*m*n;
                  };

// range: 2*min(mn,n) +3*m*n <= key
auto keyUpdateA = [=](int mn, int m, int n, int k, int l, int i)
                  {
                      return 2*mn + 3*m*n + (i-1)*m*n + m*(l-1) + (k-1);
                  };

} // namespace flens

#endif // LU_LU_TILED_MT_H
