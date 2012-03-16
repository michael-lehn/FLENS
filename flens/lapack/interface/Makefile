#CXX      = clang++
CXX      = g++

FC = gfortran

#INTEGER = long
INTEGER = int

CXXFLAGS += -Wall -Wextra
#CXXFLAGS += -Wc++0x-extensions
CXXFLAGS += -Werror
CXXFLAGS += -g
CXXFLAGS += -O3 -std=c++0x -I ../../.. -DINTEGER=$(INTEGER)
CXXFLAGS += -DWITH_REFBLAS

CXXFLAGS += -DCHECK_CXXLAPACK
#CXXFLAGS += -DUSE_NATIVE_ILAENV
#CXXFLAGS += -DNDEBUG

FCFLAGS  += -g -fno-underscoring

# trace calls to interface
CXXFLAGS += -DDEBUG_INTERFACE

CXXFLAGS += -DDEBUG_LAPACK_CALLS
CXXFLAGS += -DDEBUG_FLENS_LAPACK_CALLS


SOURCE_FILES = $(wildcard src/*.cc)

FORTRAN_FILES = $(wildcard ref_lapack/*.f) \
                $(wildcard ref_blas/*.f)

OBJECT_FILES =  $(patsubst %.cc,%.o,$(SOURCE_FILES)) \
                $(patsubst %.f,%.o,$(FORTRAN_FILES))

FLENSLAPACK = lapack_FLENS.a

all : $(OBJECT_FILES)
	ar cru $(FLENSLAPACK) $(OBJECT_FILES)

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

ref_lapack/%.o : ref_lapack/%.f
	$(FC) $(FCFLAGS) -c -o $@ $<

ref_blas/%.o : ref_blas/%.f
	$(FC) $(FCFLAGS) -c -o $@ $<

clean :
	rm -f $(OBJECT_FILES) $(FLENSLAPACK)

