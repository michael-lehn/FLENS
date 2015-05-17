FC = gfortran

FCFLAGS           += -g -fimplicit-none
FC_NOUNDERSCORING  = -fno-underscoring

SOURCE_FILES = $(filter-out ilaenv.f, $(wildcard *.f))

OBJECT_FILES_NOUNDERSCORING   =  $(filter-out xerbla_nu.o, $(patsubst %.f,%_nu.o,$(SOURCE_FILES)))

OBJECT_FILES_WITHUNDERSCORING =  $(patsubst %.f,%_wu.o,$(SOURCE_FILES))


REFLAPACK_NU = ../libreflapack_nu.a
REFLAPACK_WU = ../libreflapack_wu.a

all : $(REFLAPACK_NU) $(REFLAPACK_WU)

$(REFLAPACK_NU) : $(OBJECT_FILES_NOUNDERSCORING)
	ar crs $(REFLAPACK_NU) $(OBJECT_FILES_NOUNDERSCORING)
	ranlib $(REFLAPACK_NU)

$(REFLAPACK_WU) : $(OBJECT_FILES_WITHUNDERSCORING)
	ar crs $(REFLAPACK_WU) $(OBJECT_FILES_WITHUNDERSCORING)
	ranlib $(REFLAPACK_WU)

%_nu.o : %.f
	$(FC) $(FCFLAGS) $(FC_NOUNDERSCORING) -c -o $@ $<

%_wu.o : %.f
	$(FC) $(FCFLAGS) $(FC_WITHUNDERSCORING) -c -o $@ $<

clean :
	rm -f $(OBJECT_FILES_NOUNDERSCORING)
	rm -f $(OBJECT_FILES_WITHUNDERSCORING)
	rm -f $(REFLAPACK_NU)
	rm -f $(REFLAPACK_WU)

