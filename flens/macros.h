#ifndef FLENS_MACROS_H
#define FLENS_MACROS_H 1

#   include <cassert>

#   define ADDRESS(x) reinterpret_cast<const void *>(&x)

#   ifndef ASSERT
#       define ASSERT(x) assert(x)
#   endif //ASSERT

#endif // FLENS_MACROS_H
