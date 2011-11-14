#ifndef FLENS_LAPACK_INTERFACE_INCLUDE_CXXLAPACK_H
#define FLENS_LAPACK_INTERFACE_INCLUDE_CXXLAPACK_H 1

#include <iostream>

#ifdef LAPACK_DECL
#   undef   LAPACK_DECL
#endif 
#define  LAPACK_DECL(x)           x##_

extern "C" {
#include <flens/lapack/interface/include/lapack_test.in.h>
#include <flens/lapack/interface/include/lapack.in.h>
} // extern "C"

#undef LAPACK_DECL

#endif // FLENS_LAPACK_INTERFACE_INCLUDE_CXXLAPACK_H