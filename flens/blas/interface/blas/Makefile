CXX = g++
FC  = gfortran

CXX_SOURCE_FILES = $(wildcard *.cc) \
                   $(wildcard asum/*.cc) \
                   $(wildcard dot/*.cc) \
                   $(wildcard iamax/*.cc) \
                   $(wildcard nrm2/*.cc)
F77_SOURCE_FILES = $(wildcard asum/*.f) \
                   $(wildcard dot/*.f) \
                   $(wildcard iamax/*.f) \
                   $(wildcard nrm2/*.f)

OBJECT_FILES     = $(patsubst %.cc,%.o,$(CXX_SOURCE_FILES)) \
                   $(patsubst %.f,%.o,$(F77_SOURCE_FILES))

CXXFLAGS    += -I../../../.. -std=c++11 -Wall -Wextra -Werror
LDFLAGS     +=

FLENSBLAS  = blas_FLENS.a

all : $(OBJECT_FILES)
	ar crs $(FLENSBLAS) $(OBJECT_FILES)
	ranlib $(FLENSBLAS)

%.o : %.f
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	$(RM) $(OBJECT_FILES) $(FLENSBLAS)
