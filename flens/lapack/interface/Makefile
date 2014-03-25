CXX  = g++
#CXX  = clang++

CXXFLAGS += -g -std=c++11 -Wall -Wextra -Werror
CXXFLAGS += -I../../..

#
# FLENS-LAPACK checks its results by calling REFLAPACK and REFBALS functions
# without underscore
#
CXXFLAGS += -DWITH_REFBLAS
#CXXFLAGS += -DWITH_REFBLAS -DCXXBLAS_DEBUG

#CXXFLAGS += -DLAPACK_DEBUG

test_external_binding =

ifndef test_external_binding
#
# Call generic FLENS-LAPACK implementation and compare results with REFLAPACK
#
    CXXFLAGS += -DUSE_CXXLAPACK -DCXXLAPACK_NO_UNDERSCORE -DCHECK_CXXLAPACK \
                -DUSE_ILAENV_WITH_UNDERSCORE
else
#
# Call external LAPACK implementation (e.g. REFLAPACK).  We have define the
# macro "USE_NATIVE_ILAENV" such that the ilaenv function from the LAPACK
# test suite gets called.
#
    CXXFLAGS += -DUSE_CXXLAPACK -DCXXLAPACK_NO_UNDERSCORE \
                -DALWAYS_USE_CXXLAPACK \
                -DUSE_ILAENV_WITH_UNDERSCORE -DUSE_NATIVE_ILAENV
endif

REFLAPACK = ../../cxxlapack/netlib/libreflapack_nu.a
REFBLAS   = ../../cxxblas/netlib/libblas_REF_nu.a


SOURCE_FILES = $(wildcard src/*.cc)
#SOURCE_FILES = src/getrf.cc src/xdummy.cc

OBJECT_FILES =  $(patsubst %.cc,%.o,$(SOURCE_FILES))
FLENSLAPACK  = lapack_FLENS.a

all : $(FLENSLAPACK)

$(FLENSLAPACK) : $(OBJECT_FILES)
	ar crs $(FLENSLAPACK) $(OBJECT_FILES)
	ranlib $(FLENSLAPACK)

%.o : %.f
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean :
	rm -f $(OBJECT_FILES) $(FLENSLAPACK)
