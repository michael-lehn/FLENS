FC = gfortran
FCFLAGS = -Wall -lstdc++
BLASLIB = ../interface/libblas.a
CBLASLIB = ../interface/libcblas.a

LEVEL1TESTS = sblat1 dblat1 cblat1 zblat1

LEVEL2TESTS = sblat2 dblat2 cblat2 zblat2

LEVEL3TESTS = sblat3 dblat3 cblat3 zblat3

TARGETS = $(LEVEL1TESTS) $(LEVEL2TESTS) $(LEVEL3TESTS)


all : $(TARGETS)

$(TARGETS): % : %.o $(BLASLIB) $(CBLASLIB)
	$(FC) $(FCFLAGS) -o $@ $< $(BLASLIB) $(CBLASLIB)

%.o : %.f
	$(FC) $(FCFLAGS) -c -o $@ $<

clean :
	rm -f $(TARGETS) *.SUMM

test : $(TARGETS)
	rm -f *.SUMM *.SNAP
	./sblat1
	./dblat1
	./cblat1
	./zblat1
	./sblat2 < SBLAT2.DAT; cat SBLAT2.SUMM
	./dblat2 < DBLAT2.DAT; cat DBLAT2.SUMM
	./cblat2 < CBLAT2.DAT; cat CBLAT2.SUMM
	./zblat2 < ZBLAT2.DAT; cat ZBLAT2.SUMM
	./sblat3 < SBLAT3.DAT; cat SBLAT3.SUMM
	./dblat3 < DBLAT3.DAT; cat DBLAT3.SUMM
	./cblat3 < CBLAT3.DAT; cat CBLAT3.SUMM
	./zblat3 < ZBLAT3.DAT; cat ZBLAT3.SUMM

