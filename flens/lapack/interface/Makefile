#CXX  = g++
CXX  = clang++

CXXFLAGS += -std=c++0x -Wall -Wextra -Werror
CXXFLAGS += -I../../..

#
# FLENS-LAPACK checks its results by calling REFLAPACK and REFBALS functions
# without underscore
#
CXXFLAGS += -DWITH_REFBLAS
CXXFLAGS += -DUSE_CXXLAPACK -DCXXLAPACK_NO_UNDERSCORE -DCHECK_CXXLAPACK \
            -DUSE_ILAENV_WITH_UNDERSCORE
REFLAPACK = ../../cxxlapack/netlib/libreflapack_nu.a
REFBLAS   = ../../cxxblas/netlib/libblas_REF_nu.a


SOURCE_FILES = $(wildcard src/*.cc)
#SOURCE_FILES = src/getrf.cc src/xdummy.cc

OBJECT_FILES =  $(patsubst %.cc,%.o,$(SOURCE_FILES))
FLENSLAPACK  = lapack_FLENS.a

all : $(FLENSLAPACK)

$(FLENSLAPACK) : $(OBJECT_FILES)
	ar cru $(FLENSLAPACK) $(OBJECT_FILES)
	ranlib $(FLENSLAPACK)

%.o : %.f
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean :
	rm -f $(OBJECT_FILES) $(FLENSLAPACK)
