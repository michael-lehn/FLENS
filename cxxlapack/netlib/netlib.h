#ifndef CXXLAPACK_NETLIB_NETLIB_H
#define CXXLAPACK_NETLIB_NETLIB_H 1

#ifdef LAPACK_IMPL
#   undef   LAPACK_IMPL
#endif

#ifndef CXXLAPACK_NO_UNDERSCORE
#define     LAPACK_IMPL(x)           x##_
#else
#define     LAPACK_IMPL(x)           x
#endif

#ifndef DEBUG_CXXLAPACK
#   define  DEBUG_CXXLAPACK(x)
#endif

extern "C" {
#   include <cxxlapack/netlib/interface/lapack.in.h>
} // extern "C"

#endif //  CXXLAPACK_NETLIB_NETLIB_H
