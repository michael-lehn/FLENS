all:


check :
	$(MAKE) -C cxxblas/netlib
	$(MAKE) -C cxxlapack/netlib
	$(MAKE) -C flens/blas/interface
	$(MAKE) -C flens/lapack
	$(MAKE) -C flens/blas/interface check
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
