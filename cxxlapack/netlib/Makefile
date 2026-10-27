#
#  Reference LAPACK 3.3.1 from Netlib (used by the FLENS test suites).
#  See toolchain.mk in the FLENS root directory for configuration.
#

FLENS_ROOT ?= ../..
include $(FLENS_ROOT)/toolchain.mk

all : lapack/
	$(MAKE) -C lapack FC="$(FC) $(FC_REF)"

lapack/ :
	tar xfz lapack.tgz

clean : lapack/
	$(MAKE) -C lapack clean

distclean : clean
	$(RM) -rf lapack/

.PHONY : all clean distclean
