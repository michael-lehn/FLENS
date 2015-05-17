#ifndef LU_LU_TILED_KEYS_H
#define LU_LU_TILED_KEYS_H 1

#include <flens/examples/lu/format.h>

namespace flens {

// Some auxiliary function for generating nice task keys
auto keyLU      = [=](int i)
                  {
                      return format("LU(A_{%d,%d})", i, i);
                  };

auto keyPtA     = [=](int i, int j)
                  {
                      return format("A_{%d,%d} \\leftarrow "
                                    "P_%d^T A_{%d,%d}",
                                    i, j, i, i, j);
                  };

auto keyUpdateP = [=](int i, int bs)
                  {
                      return format("P\\big|_%d \\leftarrow P_%d + %d",
                                    i, i, (i-1)*bs);
                  };

auto keyApplyU  = [=](int i, int j)
                  {
                      return format("A_{%d,%d} \\leftarrow "
                                    "A_{%d,%d} U_{%d,%d}^{-1}",
                                    i, j, i, j, i, i);
                  };

auto keyApplyL  = [=](int i, int j)
                  {
                      return format("A_{%d,%d} \\leftarrow "
                                    "L_{%d,%d}^{-1} A_{%d,%d}",
                                    i, j, i, i, i, j);
                  };

auto keyUpdateA = [=](int i, int j, int k)
                  {
                      return format("A_{%d,%d} \\leftarrow "
                                    "A_{%d,%d} - A_{%d,%d} A_{%d,%d}",
                                    i, j, i, j, i, k, k, j);
                  };

} // namespace flens

#endif // LU_LU_TILED_MT_H
