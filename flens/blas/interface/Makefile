#
#  FLENS-BLAS with a Fortran interface (blas_FLENS.a) and the Netlib BLAS
#  test suite.
#
#    make          build blas/blas_FLENS.a
#    make check    build it and run the Netlib BLAS tests (requires FC)
#
#  See toolchain.mk in the FLENS root directory for configuration.
#

FLENS_ROOT ?= ../../..
include $(FLENS_ROOT)/toolchain.mk

FLENSBLAS = blas_FLENS.a

all : lib test/

lib :
	$(MAKE) -C blas

test/ :
	tar xfz test.tgz

#
#  The Netlib test programs always exit with status 0.  Their results are
#  collected from the output and the *.SUMM files.  The test programs are
#  relinked on every run so that they always use the current library.
#
check : lib test/
	$(MAKE) -C test clean
	$(MAKE) -C test all FC="$(FC) $(FC_COMPAT)" \
	        LDFLAGS="../blas/$(FLENSBLAS) $(CXXLIBS)"
	$(MAKE) -C test check > test/check.log 2>&1
	@cd test && \
	l1=$$(grep -c -- '----- PASS -----' check.log); \
	l23=$$(cat *.SUMM | grep -c 'PASSED THE COMPUTATIONAL'); \
	ee=$$(cat *.SUMM | grep -c 'PASSED THE TESTS OF ERROR'); \
	bad=$$(cat check.log *.SUMM | grep -ciE 'fail|suspect'); \
	echo "Netlib BLAS tests (see test/check.log, test/*.SUMM):"; \
	echo "  level 1 subprograms passed:     $$l1/40"; \
	echo "  level 2/3 computational passed: $$l23/96"; \
	echo "  error-exit tests passed:        $$ee/78"; \
	echo "  lines reporting failures:       $$bad"; \
	test $$l1 -eq 40 -a $$l23 -eq 96 -a $$ee -eq 78 -a $$bad -eq 0

clean : test/
	$(MAKE) -C blas clean
	$(MAKE) -C test clean

distclean :
	$(MAKE) -C blas clean
	$(RM) -r test/

.PHONY : all lib check clean distclean
