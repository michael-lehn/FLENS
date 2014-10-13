all:
	@echo
	@echo "FLENS is header only.  Just include -I 'PATH_TO_FLENS' in your"
	@echo "compile command, e.g.:"
	@echo
	@echo "   g++ -Wall -std=c++11 -I $(PWD) tut01-page01-example.cc"
	@echo
	@echo

install: all

flens-blas:
	$(MAKE) -C flens/blas/interface

flens-lapack:
	$(MAKE) -C flens/lapack

check-blas: flens-blas
	$(MAKE) -C flens/blas/interface check

check-lapack: flens-lapack
	$(MAKE) -C flens/lapack check

dev-check :
	$(MAKE) -C cxxblas/netlib
	$(MAKE) -C cxxlapack/netlib
	$(MAKE) -C flens/blas/interface
	$(MAKE) -C flens/lapack dev
	$(MAKE) -C flens/blas/interface dev-check
	$(MAKE) -C flens/lapack check

clean :
	$(MAKE) -C cxxblas/netlib clean
	$(MAKE) -C cxxlapack/netlib clean
	$(MAKE) -C flens/lapack clean

distclean :
	$(MAKE) -C cxxblas/netlib distclean
	$(MAKE) -C cxxlapack/netlib distclean
	$(MAKE) -C flens/blas/interface distclean
	$(MAKE) -C flens/lapack distclean
