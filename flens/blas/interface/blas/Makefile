CXX = g++

SOURCE_FILES = $(wildcard *.cc)
OBJECT_FILES = $(patsubst %.cc,%.o,$(SOURCE_FILES))
TARGETS      = $(patsubst %.cc,%,$(SOURCE_FILES))

CXXFLAGS    = -I../../../.. -std=c++11
LDFLAGS     =

FLENSBLAS  = blas_FLENS.a

all : $(OBJECT_FILES)
	ar crs $(FLENSBLAS) $(OBJECT_FILES)
	ranlib $(FLENSBLAS)

clean:
	$(RM) $(TARGETS) $(OBJECT_FILES)
