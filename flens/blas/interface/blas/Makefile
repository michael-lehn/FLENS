CXX = g++-4.8
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


CXXFLAGS    += -I../../../.. -std=c++11 -Wall -Wextra -Werror
#CXXFLAGS    += -DBLAS_DEBUG
LDFLAGS     +=

ifndef NO_UNDERSCORE
    OBJECT_FILES     = $(patsubst %.cc,%.o,$(CXX_SOURCE_FILES)) \
                       $(patsubst %.f,%.o,$(F77_SOURCE_FILES))
    FLENSBLAS  = blas_FLENS.a
else
    OBJECT_FILES     = $(patsubst %.cc,%_nu.o,$(CXX_SOURCE_FILES)) \
                       $(patsubst %.f,%_nu.o,$(F77_SOURCE_FILES))
    FCFLAGS         += -fno-underscoring
    CXXFLAGS        += -DNOCHANGE
    FLENSBLAS        = blas_FLENS_nu.a
endif

all : $(FLENSBLAS)

$(FLENSBLAS) : $(OBJECT_FILES)
	ar crs $(FLENSBLAS) $(OBJECT_FILES)
	ranlib $(FLENSBLAS)

%.o : %.f
	$(FC) $(FCFLAGS) -c -o $@ $<

%_nu.o : %.f
	$(FC) $(FCFLAGS) -c -o $@ $<

%.o : %.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c -o $@ $<

%_nu.o : %.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c -o $@ $<

clean:
	$(RM) $(OBJECT_FILES) $(FLENSBLAS)
