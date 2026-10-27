#
#  Toolchain configuration shared by the FLENS Makefiles.
#
#  All variables can be set on the command line or in the environment, e.g.
#
#    make check CXX=clang++ FC=gfortran-14 CXXSTD=c++20
#    make check CXXFLAGS="-O0 -g" WERROR=1
#
#  CXX        C++ compiler (default: c++)
#  CXXSTD     C++ standard (default: c++11, the minimum required by FLENS)
#  CXXFLAGS   C++ compiler flags (default: -O2 -Wall -Wextra)
#  WERROR     if non-empty, add -Werror to CXXFLAGS
#  FC         Fortran compiler (default: gfortran), needed for the Netlib
#             reference implementations and test suites
#  FC_COMPAT  extra Fortran flags required by the old Netlib sources
#             (default: -fallow-argument-mismatch if FC accepts it; gfortran
#             10 and later reject argument mismatches otherwise)
#  CXXLIBS    C++ runtime library for linking C++ objects with FC
#             (default: -lc++ for clang on macOS, -lstdc++ otherwise)
#  FC_REF     extra Fortran flags for the Netlib reference libraries
#             (default: FC_COMPAT and -ffp-contract=off if FC accepts it)
#  DEV_CXXFLAGS
#             extra C++ flags for the development build of FLENS-LAPACK
#             used by 'make dev-check', which compares results bit by bit
#             with the reference implementation (default: -ffp-contract=off
#             and -fcx-fortran-rules if CXX accepts them)
#

ifeq ($(origin CXX),default)
    CXX = c++
endif

ifeq ($(origin FC),default)
    FC = gfortran
endif

CXXSTD   ?= c++11
CXXFLAGS ?= -O2 -Wall -Wextra

ifneq ($(strip $(WERROR)),)
    CXXFLAGS += -Werror
endif

ifeq ($(origin FC_COMPAT),undefined)
    FC_COMPAT := $(shell $(FC) -fallow-argument-mismatch -x f77 -c /dev/null \
                         -o /dev/null >/dev/null 2>&1 \
                         && echo -fallow-argument-mismatch)
endif

ifeq ($(origin FC_REF),undefined)
    FC_REF := $(FC_COMPAT) \
              $(shell $(FC) -ffp-contract=off -x f77 -c /dev/null \
                      -o /dev/null >/dev/null 2>&1 && echo -ffp-contract=off)
endif

ifeq ($(origin DEV_CXXFLAGS),undefined)
    DEV_CXXFLAGS := \
        $(shell $(CXX) -ffp-contract=off -x c++ -c /dev/null \
                -o /dev/null >/dev/null 2>&1 && echo -ffp-contract=off) \
        $(shell $(CXX) -fcx-fortran-rules -x c++ -c /dev/null \
                -o /dev/null >/dev/null 2>&1 && echo -fcx-fortran-rules)
endif

ifeq ($(origin CXXLIBS),undefined)
    ifneq ($(findstring clang,$(shell $(CXX) --version 2>/dev/null)),)
        ifeq ($(shell uname -s),Darwin)
            LIBCXX_DIR := $(wildcard $(shell $(CXX) -print-resource-dir \
                                     2>/dev/null)/../../c++)
            CXXLIBS := -lc++
            ifneq ($(LIBCXX_DIR),)
                CXXLIBS := -L$(LIBCXX_DIR) -Wl,-rpath,$(LIBCXX_DIR) -lc++
            endif
        else
            CXXLIBS := -lstdc++
        endif
    else
        CXXLIBS := -lstdc++
    endif
endif

export CXX FC CXXSTD CXXFLAGS FC_COMPAT FC_REF DEV_CXXFLAGS CXXLIBS
