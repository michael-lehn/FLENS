CXX      = clang++
#CXX      = g++

CXXFLAGS = -I ../.. -Wall

TEST_ROW_MAJOR = U
#TEST_ROW_MAJOR = D


CBLAS_INT = int
BU        = _

BLAS_iSDCZ  = amax
BLAS_sCdZ   = nrm2 asum
BLAS_SD     = dot rot rotg nrm2 asum \
              symv sbmv spmv syr spr syr2 spr2
BLAS_SDCZ   = copy swap axpy scal \
              gemv gbmv trmv tbmv tpmv trsv tbsv tpsv ger gemm trmm trsm \
              symm syrk syr2k
BLAS_CZ     = dotu dotc gerc geru hemv hbmv hpmv her hpr her2 hpr2 hemm \
              herk her2k
BLAS_CsZd   = scal


BLAS_OBJ =  $(patsubst %,is%.o,$(BLAS_iSDCZ)) \
            $(patsubst %,id%.o,$(BLAS_iSDCZ)) \
            $(patsubst %,ic%.o,$(BLAS_iSDCZ)) \
            $(patsubst %,iz%.o,$(BLAS_iSDCZ)) \
            $(patsubst %,sc%.o,$(BLAS_sCdZ)) \
            $(patsubst %,dz%.o,$(BLAS_sCdZ)) \
            $(patsubst %,s%.o,$(BLAS_SDCZ) $(BLAS_SD)) \
            $(patsubst %,d%.o,$(BLAS_SDCZ) $(BLAS_SD)) \
            $(patsubst %,c%.o,$(BLAS_SDCZ) $(BLAS_CZ)) \
            $(patsubst %,z%.o,$(BLAS_SDCZ) $(BLAS_CZ)) \
            $(patsubst %,cs%.o,$(BLAS_CsZd)) \
            $(patsubst %,zd%.o,$(BLAS_CsZd))

CBLAS_OBJ = $(patsubst %,cblas_%,$(BLAS_OBJ))

BLAS_OPTION = -DBLAS_NAME=<tchar>$(*F)$(BU) -DCBLAS_NAME=cblas_<tchar>$(*F) \
              -<complex1>COMPLEX_FLOAT1 -<complex2>COMPLEX_FLOAT2 \
              -DCBLAS_FLOAT=<type> -DCBLAS_INT=$(CBLAS_INT) \
              -DBLAS_NAME_STR="\"<tchar>$(*F)      \"" \
              -$(TEST_ROW_MAJOR)TEST_ROW_MAJOR

is_BLAS_OPTION = $(subst \
                    <type>,float,$(subst \
                    <tchar>,is,$(subst \
                    <complex1>,U,$(subst \
                    <complex2>,U,$(BLAS_OPTION)))))

id_BLAS_OPTION = $(subst \
                    <type>,double,$(subst \
                    <tchar>,id,$(subst \
                    <complex1>,U,$(subst \
                    <complex2>,U,$(BLAS_OPTION)))))

ic_BLAS_OPTION = $(subst \
                    <type>,float,$(subst \
                    <tchar>,ic,$(subst \
                    <complex1>,D,$(subst \
                    <complex2>,D,$(BLAS_OPTION)))))

iz_BLAS_OPTION = $(subst \
                    <type>,double,$(subst \
                    <tchar>,iz,$(subst \
                    <complex1>,D,$(subst \
                    <complex2>,D,$(BLAS_OPTION)))))


sc_BLAS_OPTION = $(subst \
                    <type>,float,$(subst \
                    <tchar>,sc,$(subst \
                    <complex1>,D,$(subst \
                    <complex2>,U,$(BLAS_OPTION)))))

dz_BLAS_OPTION = $(subst \
                    <type>,double,$(subst \
                    <tchar>,dz,$(subst \
                    <complex1>,D,$(subst \
                    <complex2>,U,$(BLAS_OPTION)))))


cs_BLAS_OPTION = $(subst \
                    <type>,float,$(subst \
                    <tchar>,cs,$(subst \
                    <complex1>,D,$(subst \
                    <complex2>,U,$(BLAS_OPTION)))))

zd_BLAS_OPTION = $(subst \
                    <type>,double,$(subst \
                    <tchar>,zd,$(subst \
                    <complex1>,D,$(subst \
                    <complex2>,U,$(BLAS_OPTION)))))


s_BLAS_OPTION = $(subst \
                    <type>,float,$(subst \
                    <tchar>,s,$(subst \
                    <complex1>,U,$(subst \
                    <complex2>,U,$(BLAS_OPTION)))))

d_BLAS_OPTION = $(subst \
                    <type>,double,$(subst \
                    <tchar>,d,$(subst \
                    <complex1>,U,$(subst \
                    <complex2>,U,$(BLAS_OPTION)))))

c_BLAS_OPTION = $(subst \
                    <type>,float,$(subst \
                    <tchar>,c,$(subst \
                    <complex1>,D,$(subst \
                    <complex2>,D,$(BLAS_OPTION)))))

z_BLAS_OPTION = $(subst \
                    <type>,double,$(subst \
                    <tchar>,z,$(subst \
                    <complex1>,D,$(subst \
                    <complex2>,D,$(BLAS_OPTION)))))


all : $(BLAS_OBJ) $(CBLAS_OBJ)
	ar ru libblas.a $(BLAS_OBJ)
	ar ru libcblas.a $(CBLAS_OBJ)

is%.o : %.cc
	$(CXX) $(CXXFLAGS) $(is_BLAS_OPTION) -DCREATE_BLAS -c -o $@ $<

id%.o : %.cc
	$(CXX) $(CXXFLAGS) $(id_BLAS_OPTION) -DCREATE_BLAS -c -o $@ $<

ic%.o : %.cc
	$(CXX) $(CXXFLAGS) $(ic_BLAS_OPTION) -DCREATE_BLAS -c -o $@ $<

iz%.o : %.cc
	$(CXX) $(CXXFLAGS) $(iz_BLAS_OPTION) -DCREATE_BLAS -c -o $@ $<


sc%.o : %.cc
	$(CXX) $(CXXFLAGS) $(sc_BLAS_OPTION) -DCREATE_BLAS -c -o $@ $<

dz%.o : %.cc
	$(CXX) $(CXXFLAGS) $(dz_BLAS_OPTION) -DCREATE_BLAS -c -o $@ $<


cs%.o : %.cc
	$(CXX) $(CXXFLAGS) $(cs_BLAS_OPTION) -DCREATE_BLAS -c -o $@ $<

zd%.o : %.cc
	$(CXX) $(CXXFLAGS) $(zd_BLAS_OPTION) -DCREATE_BLAS -c -o $@ $<


s%.o : %.cc
	$(CXX) $(CXXFLAGS) $(s_BLAS_OPTION) -DCREATE_BLAS -c -o $@ $<

d%.o : %.cc
	$(CXX) $(CXXFLAGS) $(d_BLAS_OPTION) -DCREATE_BLAS -c -o $@ $<

c%.o : %.cc
	$(CXX) $(CXXFLAGS) $(c_BLAS_OPTION) -DCREATE_BLAS -c -o $@ $<

z%.o : %.cc
	$(CXX) $(CXXFLAGS) $(z_BLAS_OPTION) -DCREATE_BLAS -c -o $@ $<


cblas_is%.o : %.cc
	$(CXX) $(CXXFLAGS) $(is_BLAS_OPTION) -DCREATE_CBLAS -c -o $@ $<

cblas_id%.o : %.cc
	$(CXX) $(CXXFLAGS) $(id_BLAS_OPTION) -DCREATE_CBLAS -c -o $@ $<

cblas_ic%.o : %.cc
	$(CXX) $(CXXFLAGS) $(ic_BLAS_OPTION) -DCREATE_CBLAS -c -o $@ $<

cblas_iz%.o : %.cc
	$(CXX) $(CXXFLAGS) $(iz_BLAS_OPTION) -DCREATE_CBLAS -c -o $@ $<


cblas_sc%.o : %.cc
	$(CXX) $(CXXFLAGS) $(sc_BLAS_OPTION) -DCREATE_CBLAS -c -o $@ $<

cblas_dz%.o : %.cc
	$(CXX) $(CXXFLAGS) $(dz_BLAS_OPTION) -DCREATE_CBLAS -c -o $@ $<


cblas_cs%.o : %.cc
	$(CXX) $(CXXFLAGS) $(cs_BLAS_OPTION) -DCREATE_CBLAS -c -o $@ $<

cblas_zd%.o : %.cc
	$(CXX) $(CXXFLAGS) $(zd_BLAS_OPTION) -DCREATE_CBLAS -c -o $@ $<


cblas_s%.o : %.cc
	$(CXX) $(CXXFLAGS) $(s_BLAS_OPTION) -DCREATE_CBLAS -c -o $@ $<

cblas_d%.o : %.cc
	$(CXX) $(CXXFLAGS) $(d_BLAS_OPTION) -DCREATE_CBLAS -c -o $@ $<

cblas_c%.o : %.cc
	$(CXX) $(CXXFLAGS) $(c_BLAS_OPTION) -DCREATE_CBLAS -c -o $@ $<

cblas_z%.o : %.cc
	$(CXX) $(CXXFLAGS) $(z_BLAS_OPTION) -DCREATE_CBLAS -c -o $@ $<

clean :
	rm -f *.o *.a

