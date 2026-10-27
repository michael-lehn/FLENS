FLENS_ROOT ?= ../../../..
include $(FLENS_ROOT)/toolchain.mk

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

ALL_CXXFLAGS = -std=$(CXXSTD) $(CXXFLAGS) -I$(FLENS_ROOT) -MMD -MP

FLENSBLAS  = blas_FLENS.a

all : $(FLENSBLAS)

$(FLENSBLAS) : $(OBJECT_FILES)
	$(RM) $(FLENSBLAS)
	ar crs $(FLENSBLAS) $(OBJECT_FILES)
	ranlib $(FLENSBLAS)

%.o : %.cc
	$(CXX) $(ALL_CXXFLAGS) -c -o $@ $<

%.o : %.f
	$(FC) $(FCFLAGS) -c -o $@ $<

#
# Check that the library exports exactly the routines listed in exports.txt
#
check-exports : $(FLENSBLAS)
	$(FLENS_ROOT)/check-exports.sh $(FLENSBLAS) exports.txt

clean:
	$(RM) $(OBJECT_FILES) $(OBJECT_FILES:.o=.d) $(FLENSBLAS)

.PHONY : all check-exports clean

-include $(OBJECT_FILES:.o=.d)
