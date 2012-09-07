all : blas/blas_FLENS.a test/

blas/blas_FLENS.a :
	$(MAKE) -C blas

test/ :
	tar xfz test.tgz

check : test/
	$(MAKE) -C test
	$(MAKE) -C test check

clean : test/
	$(MAKE) -C blas clean
	$(MAKE) -C test clean

distclean :
	$(MAKE) -C blas clean
	$(RM) -rf test/
