DLINTST := $(TOPDIR)/xdlintst
TMGLIB  := $(TOPDIR)/libtmglib.a

CXX        = g++-4.8
FC         =  gfortran
FFLAGS     += -fimplicit-none -O3



FCLIBSEARCH := $(shell $(FC) -print-search-dirs  \
                      | grep libraries          \
                      | sed 's,libraries:.*=,,' \
                      | sed 's,:, ,g')
LDF77FLAGS  += $(shell for path in $(FCLIBSEARCH); \
                      do                            \
                          if [ -d "$$path" ]; then    \
                              echo "-L $$path "; \
                          fi;                       \
                      done)
LDF77LIBS   += -lgfortran
