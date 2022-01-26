CXX  = clang++

CXXFLAGS += -g -std=c++11 -Wall -Wextra -Werror
CXXFLAGS += -I../../..

CXXFLAGS += -DWITH_REFBLAS
#CXXFLAGS += -DTRACEBACK_ASSERT

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
	ar crs $(FLENSLAPACK) $(OBJECT_FILES)
	ranlib $(FLENSLAPACK)

$(FLENSLAPACK_DEV) : $(OBJECT_FILES_DEV)
	ar crs $(FLENSLAPACK) $(OBJECT_FILES)
	ranlib $(FLENSLAPACK)

%.o : %.cc
	$(CXX) $(CXXFLAGS) $(CHECK_CXXFLAGS) -c -o $@ $<

%.dev_o : %.cc
	$(CXX) $(CXXFLAGS) $(DEV_CHECK_CXXFLAGS) -c -o $@ $<

clean :
	rm -f $(OBJECT_FILES) $(OBJECT_FILES_DEV) $(FLENSLAPACK) $(FLENSLAPACK_DEV)
