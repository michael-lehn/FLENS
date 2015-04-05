FC = gfortran

SOURCE_FILES = $(wildcard *.f)
OBJECT_FILES = $(patsubst %.f,%.o,$(SOURCE_FILES))
TARGETS      = $(patsubst %.f,%,$(SOURCE_FILES))

FLENSBLAS = ../blas/blas_FLENS.a

#LDFLAGS     += $(FLENSBLAS) -L. -lrefblas -lstdc++
LDFLAGS     += $(FLENSBLAS) -L. -lstdc++

all : $(TARGETS)

%: %.f
	$(FC) $< $(LDFLAGS) -o $@
check:
	$(RM) *.SUMM
	./sblat1
	./dblat1
	./cblat1
	./zblat1
	./sblat2 < sblat2.dat
	./dblat2 < dblat2.dat
	./cblat2 < cblat2.dat
	./zblat2 < zblat2.dat
	./sblat3 < sblat3.dat
	./dblat3 < dblat3.dat
	./cblat3 < cblat3.dat
	./zblat3 < zblat3.dat

clean:
	$(RM) $(TARGETS) $(OBJECT_FILES) *.SUMM
