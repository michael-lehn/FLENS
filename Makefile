#
#  FLENS is header only.  This Makefile builds and runs the tests:
#
#    make check                  regression tests, Netlib BLAS and LAPACK
#                                test suites (the latter require a Fortran
#                                compiler)
#    make check-regression       regression tests only (C++ only)
#    make check-blas             Netlib BLAS tests for FLENS-BLAS
#    make check-lapack           Netlib LAPACK tests for FLENS-LAPACK
#    make check-lapack-extended  all double/double complex LAPACK tests
#    make dev-check              BLAS tests and LAPACK tests with FLENS-LAPACK
#                                routines comparing their results with the
#                                reference implementation
#
#  check-blas (blas_FLENS.a) is checked first: the LAPACK targets require
#  it, and with 'make check' the regression tests run only after it passed.
#  A standalone 'make check-regression' does not need a Fortran compiler and
#  therefore does not require check-blas.
#
#  Compilers and flags are configured in toolchain.mk, e.g.
#
#    make check CXX=clang++ CXXSTD=c++20 FC=gfortran
#

include toolchain.mk

all:
	@echo
	@echo "FLENS is header only.  Just include -I 'PATH_TO_FLENS' in your"
	@echo "compile command, e.g.:"
	@echo
	@echo "   g++ -Wall -std=c++11 -I $(CURDIR) tut01-page01-example.cc"
	@echo
	@echo "Run the tests with 'make check' (see Makefile for details)."
	@echo

install: all

flens-blas:
	$(MAKE) -C flens/blas/interface

flens-lapack:
	$(MAKE) -C flens/lapack

check: check-blas check-regression check-lapack

check-regression: | $(if $(filter check,$(MAKECMDGOALS)),check-blas)
	$(MAKE) -C flens/test/regression check

check-blas:
	$(MAKE) -C flens/blas/interface check

check-lapack: check-blas
	$(MAKE) -C flens/lapack check

#
#  The LAPACK targets share the directory flens/lapack/test.  If several of
#  them are requested (e.g. make -j check dev-check) they must not run at the
#  same time.
#
LAPACK_GOALS = $(filter check check-lapack,$(MAKECMDGOALS))

check-lapack-extended: check-blas | $(LAPACK_GOALS)
	$(MAKE) -C flens/lapack check-extended

dev-check: check-blas | $(LAPACK_GOALS) \
                       $(filter check-lapack-extended,$(MAKECMDGOALS))
	$(MAKE) -C flens/lapack dev-check

clean :
	$(MAKE) -C cxxblas/netlib clean
	$(MAKE) -C cxxlapack/netlib clean
	$(MAKE) -C flens/blas/interface clean
	$(MAKE) -C flens/lapack clean
	$(MAKE) -C flens/test/regression clean

distclean :
	$(MAKE) -C cxxblas/netlib distclean
	$(MAKE) -C cxxlapack/netlib distclean
	$(MAKE) -C flens/blas/interface distclean
	$(MAKE) -C flens/lapack distclean
	$(MAKE) -C flens/test/regression clean

.PHONY: all install flens-blas flens-lapack check check-regression \
        check-blas check-lapack check-lapack-extended dev-check clean distclean
