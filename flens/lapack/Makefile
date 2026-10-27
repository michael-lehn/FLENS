#
#  FLENS-LAPACK with a Fortran interface (interface/lapack_FLENS.a) and the
#  Netlib LAPACK 3.3.1 test suite (double and double complex precision).
#
#    make                 build interface/lapack_FLENS.a
#    make check           run the LAPACK tests zsg, dtest, ztest, dnep, znep
#    make check-extended  run all double and double complex LAPACK tests
#    make dev-check       like check, but with interface/lapack_FLENS_dev.a,
#                         in which each FLENS-LAPACK routine compares its
#                         results with the reference implementation
#
#  The results are evaluated by check-results.sh; known deviations are
#  listed in check-expected.txt.  See toolchain.mk in the FLENS root
#  directory for configuration.
#

FLENS_ROOT ?= ../..
include $(FLENS_ROOT)/toolchain.mk

CHECK_TESTS    = zsg.out dtest.out ztest.out dnep.out znep.out

EXTENDED_TESTS = dtest.out ztest.out \
                 dbak.out dbal.out dbb.out dcsd.out dec.out ded.out \
                 dgbak.out dgbal.out dgd.out dgg.out dglm.out dgqr.out \
                 dgsv.out dlse.out dnep.out dsb.out dsep.out dsg.out \
                 dsvd.out \
                 zbak.out zbal.out zbb.out zcsd.out zec.out zed.out \
                 zgbak.out zgbal.out zgd.out zgg.out zglm.out zgqr.out \
                 zgsv.out zlse.out znep.out zsb.out zsep.out zsg.out \
                 zsvd.out

TEST_PROGRAMS  = xlintstd xlintstz xeigtstd xeigtstz

TEST_VARS      = "FORTRAN=$(FC) -fimplicit-none -g $(FC_COMPAT)" \
                 "LOADER=$(FC) -g $(CXXLIBS)" \
                 "NOOPT=-g -O0 $(FC_COMPAT)"

all : lib

lib :
	$(MAKE) -C interface

dev :
	$(MAKE) -C interface dev

netlib :
	$(MAKE) -C $(FLENS_ROOT)/cxxblas/netlib
	$(MAKE) -C $(FLENS_ROOT)/cxxlapack/netlib

test/ :
	tar xfz test.tgz

#
#  run-tests: remove old test programs and results (the test programs do not
#  depend on the library and would not be relinked otherwise), build the
#  test programs one after another (they share object files), then run the
#  tests (in parallel with make -j) and evaluate the results.
#
define run-tests
	$(MAKE) -C test cleanup
	$(MAKE) -C test/MATGEN $(TEST_VARS)
	$(MAKE) -C test -j1 $(TEST_PROGRAMS) $(TEST_VARS) $(1)
	-$(MAKE) -C test -k $(2) $(TEST_VARS) $(1)
	./check-results.sh check-expected.txt $(addprefix test/,$(2))
endef

check : lib netlib test/
	$(call run-tests,,$(CHECK_TESTS))

check-extended : lib netlib test/
	$(call run-tests,,$(EXTENDED_TESTS))

dev-check : dev netlib test/
	$(call run-tests,LAPACKLIB=interface/lapack_FLENS_dev.a,$(CHECK_TESTS))

clean : test/
	$(MAKE) -C interface clean
	$(MAKE) -C test cleanup
	$(RM) interface/tmglib_FLENS.a

distclean :
	$(MAKE) -C interface clean
	$(RM) interface/tmglib_FLENS.a
	$(RM) -rf test/

.PHONY : all lib dev netlib check check-extended dev-check clean distclean
