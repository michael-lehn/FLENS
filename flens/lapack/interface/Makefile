CXX      = clang++
#CXX      = icc

FC = gfortran

INT = int

CXXFLAGS += -DWITH_REFBLAS -Werror

CXXFLAGS += -O3 -std=c++0x -I ../../.. -Wall -DINT=$(INT) -DDOUBLE -Wc++0x-extensions
FCFLAGS  += -fno-underscoring

# trace calls to interface
#CXXFLAGS += -DDEBUG_INTERFACE


CXXFLAGS += -DCMP_W_GESV_REF
CXXFLAGS += -DCMP_W_GETRF_REF
CXXFLAGS += -DCMP_W_GETF2_REF
CXXFLAGS += -DCMP_W_GEQRF_REF
CXXFLAGS += -DCMP_W_ORMQR_REF

CXXFLAGS += -DCMP_W_ORGQR_REF

CXXFLAGS += -DUSE_LASWP_REF


CXXFLAGS += -DUSE_LABAD_REF
CXXFLAGS += -DUSE_LACPY_REF
CXXFLAGS += -DUSE_LANGE_REF
CXXFLAGS += -DUSE_LANTR_REF
CXXFLAGS += -DUSE_LAQGE_REF
CXXFLAGS += -DUSE_LARNV_REF
CXXFLAGS += -DUSE_LARTG_REF
CXXFLAGS += -DUSE_LASET_REF
CXXFLAGS += -DUSE_LASSQ_REF
CXXFLAGS += -DUSE_LATRS_REF



CXXFLAGS += -DUSE_GECON_REF
CXXFLAGS += -DUSE_GEEQU_REF
CXXFLAGS += -DUSE_GERFS_REF
CXXFLAGS += -DUSE_GESVX_REF
CXXFLAGS += -DUSE_GETRI_REF




CXXFLAGS += -DUSE_GEQRFP_REF
CXXFLAGS += -DUSE_GEQR2_REF
CXXFLAGS += -DUSE_GEQR2P_REF
CXXFLAGS += -DUSE_ORG2R_REF
CXXFLAGS += -DUSE_ORM2R_REF

CXXFLAGS += -DUSE_LANSY_REF



SOURCE_FILES = $(wildcard *.cc) ../../../cxxblas/drivers/drivers.cc 

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

