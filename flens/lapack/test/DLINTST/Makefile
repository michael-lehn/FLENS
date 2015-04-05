TOPDIR = ..

include $(TOPDIR)/make.inc

SOURCE_FILES = $(wildcard *.f) $(wildcard ../COMMON/*.f) $(wildcard ../NOTPORTED/*.f)
OBJECT_FILES = $(patsubst %.f, %.o, $(SOURCE_FILES))

LDLIBS      +=  -lstdc++

LAPACK_LIB       = ../../interface/lapack_FLENS.a
NETLIB_LIB       = ../../../../cxxlapack/netlib/libreflapack_nu.a
FLENSBLAS_NU_LIB = ../../../blas/interface/blas/blas_FLENS_nu.a
F77BLAS_LIB      = ../../../blas/interface/blas/blas_FLENS.a



all: $(DLINTST)

$(DLINTST) : $(OBJECT_FILES) $(TMGLIB) $(LAPACK_LIB) $(F77BLAS_LIB) $(NETLIB_LIB) $(FLENSBLAS_NU_LIB)
	$(CXX) $(LDFLAGS) $(LDF77FLAGS) $(LDF77LIBS) -o $@ $^ $(LDLIBS)

clean:
	$(RM) $(OBJECT_FILES) $(DLINTST)
