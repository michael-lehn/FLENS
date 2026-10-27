FLENS_ROOT ?= ../../..
include $(FLENS_ROOT)/toolchain.mk

ALL_CXXFLAGS = -std=$(CXXSTD) $(CXXFLAGS) -g -I$(FLENS_ROOT) -MMD -MP

ALL_CXXFLAGS += -DWITH_REFBLAS
#ALL_CXXFLAGS += -DTRACEBACK_ASSERT

#
# Call generic FLENS-LAPACK implementation and compare results with REFLAPACK
#
DEV_CHECK_CXXFLAGS += -DUSE_CXXLAPACK -DCXXLAPACK_NO_UNDERSCORE -DCHECK_CXXLAPACK \
                      -DUSE_ILAENV_WITH_UNDERSCORE

#
# Call external LAPACK implementation (e.g. REFLAPACK).  We have define the
# macro "USE_NATIVE_ILAENV" such that the ilaenv function from the LAPACK
# test suite gets called.
#
CHECK_CXXFLAGS += -DUSE_CXXLAPACK -DCXXLAPACK_NO_UNDERSCORE \
                  -DUSE_ILAENV_WITH_UNDERSCORE -DUSE_NATIVE_ILAENV

REFLAPACK = ../../cxxlapack/netlib/libreflapack_nu.a
REFBLAS   = ../../cxxblas/netlib/libblas_REF_nu.a


SOURCE_FILES = $(wildcard src/*.cc)
#SOURCE_FILES = src/getrf.cc src/xdummy.cc

OBJECT_FILES     =  $(patsubst %.cc,%.o,$(SOURCE_FILES))
OBJECT_FILES_DEV =  $(patsubst %.cc,%.dev_o,$(SOURCE_FILES))
FLENSLAPACK      = lapack_FLENS.a
FLENSLAPACK_DEV  = lapack_FLENS_dev.a

all : $(FLENSLAPACK)

dev : $(FLENSLAPACK_DEV)

$(FLENSLAPACK) : $(OBJECT_FILES)
	$(RM) $(FLENSLAPACK)
	ar crs $(FLENSLAPACK) $(OBJECT_FILES)
	ranlib $(FLENSLAPACK)

$(FLENSLAPACK_DEV) : $(OBJECT_FILES_DEV)
	$(RM) $(FLENSLAPACK_DEV)
	ar crs $(FLENSLAPACK_DEV) $(OBJECT_FILES_DEV)
	ranlib $(FLENSLAPACK_DEV)

%.o : %.cc
	$(CXX) $(ALL_CXXFLAGS) $(CHECK_CXXFLAGS) -c -o $@ $<

%.dev_o : %.cc
	$(CXX) $(ALL_CXXFLAGS) $(DEV_CXXFLAGS) $(DEV_CHECK_CXXFLAGS) \
	       -MF $*.dev_d -MT $@ -c -o $@ $<

#
# Check that the libraries export exactly the routines listed in exports.txt
#
check-exports : $(FLENSLAPACK)
	$(FLENS_ROOT)/check-exports.sh $(FLENSLAPACK) exports.txt

check-exports-dev : $(FLENSLAPACK_DEV)
	$(FLENS_ROOT)/check-exports.sh $(FLENSLAPACK_DEV) exports.txt

clean :
	rm -f $(OBJECT_FILES) $(OBJECT_FILES_DEV) $(FLENSLAPACK) $(FLENSLAPACK_DEV)
	rm -f $(OBJECT_FILES:.o=.d) $(OBJECT_FILES_DEV:.dev_o=.dev_d)

.PHONY : all dev check-exports check-exports-dev clean

-include $(OBJECT_FILES:.o=.d) $(OBJECT_FILES_DEV:.dev_o=.dev_d)
