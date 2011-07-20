CXX      = clang++
#CXX      = g++

FC = gfortran

INT = int

# used in gesv.cc
CXXFLAGS += -DUSE_GECON_REF
CXXFLAGS += -DUSE_GEEQU_REF
CXXFLAGS += -DUSE_GERFS_REF
#CXXFLAGS += -DUSE_GESV_REF
CXXFLAGS += -DUSE_GESVX_REF
#CXXFLAGS += -DUSE_GETF2_REF
#CXXFLAGS += -DUSE_GETRF_REF
CXXFLAGS += -DUSE_GETRI_REF
#CXXFLAGS += -DUSE_GETRS_REF

# used in utils.cc
CXXFLAGS += -DUSE_LAMCH_REF

# used in la.cc
CXXFLAGS += -DUSE_LABAD_REF
CXXFLAGS += -DUSE_LACPY_REF
CXXFLAGS += -DUSE_LANGE_REF
CXXFLAGS += -DUSE_LANTR_REF
CXXFLAGS += -DUSE_LARNV_REF
CXXFLAGS += -DUSE_LARTG_REF
CXXFLAGS += -DUSE_LASET_REF
#CXXFLAGS += -DUSE_LASWP_REF
CXXFLAGS += -DUSE_LAQGE_REF




CXXFLAGS += -I ../.. -Wall -DINT=$(INT) -DSINGLE

FCFLAGS += -fno-underscoring

TEST_ROW_MAJOR = U
#TEST_ROW_MAJOR = D


BU  = _

SOURCE_FILES = utils.cc gesv.cc la.cc

FORTRAN_FILES = $(wildcard ref/*.f)

OBJECT_FILES =  $(patsubst %.cc,%.o,$(SOURCE_FILES)) \
                $(patsubst %.f,%.o,$(FORTRAN_FILES))


all : $(OBJECT_FILES)
	ar ru lapack_CXX.a $(OBJECT_FILES)
	mv lapack_CXX.a ..

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

ref/%.o : ref/%.f
	$(FC) $(FCFLAGS) -c -o $@ $<

clean :
	rm -f *.o *.a

