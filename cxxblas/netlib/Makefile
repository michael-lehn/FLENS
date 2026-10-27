#
#  Reference BLAS and CBLAS from Netlib (used by the FLENS test suites).
#  See toolchain.mk in the FLENS root directory for configuration.
#

FLENS_ROOT ?= ../..
include $(FLENS_ROOT)/toolchain.mk

all : blas/ cblas/
	$(MAKE) -C blas FC="$(FC) $(FC_REF)"
	$(MAKE) -C cblas

blas/ :
	tar xfz blas.tgz

cblas/ :
	tar xfz cblas.tgz

clean : blas/ cblas/
	$(MAKE) -C blas clean
	$(MAKE) -C cblas clean

distclean : clean
	$(RM) -rf blas/
	$(RM) -rf cblas/

.PHONY : all clean distclean
