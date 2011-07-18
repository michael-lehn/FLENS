ifndef FLENS_HOME
FLENS_HOME = $(PWD)/..
endif

include ../Makefile.common

DYLIB   = libcxxblas.$(DYLIB_EXT)
OBJDIR  = .obj/
DEPDIR  = .dep/
SOURCES := $(wildcard *.cc)
HEADERS := $(wildcard *.h) $(wildcard *.tcc)
OBJECTS = $(addprefix $(OBJDIR),$(SOURCES:.cc=.o))
DEPS    = $(addprefix $(DEPDIR),$(SOURCES:.cc=.d))

.PHONY : all
all: $(OBJDIR) $(DEPDIR) $(OBJECTS)
	$(CXX) $(CXXDYLIB) $(CXXFLAGS) -o $(DYLIB) $(OBJDIR)*.o $(LDFLAGS)
	cp $(DYLIB) $(FLENS_HOME)
	
ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif

.PHONY : install
install: $(OBJECTS)
	$(INSTALL) $(HEADERS) $(INCDIR)
	
$(OBJDIR):
	$(MKDIR) $(OBJDIR)
	
$(DEPDIR):
	$(MKDIR) $(DEPDIR)
	
$(OBJDIR)%.o: %.cc
	$(CXX) $(INCDIRS) $(CXXFLAGS) -o $@ -c $<

$(DEPDIR)%.d: %.cc $(DEPDIR)
	$(CXX) -MM $(INCDIRS) $(CXXFLAGS) $< \
	| sed 's/$*.o/$(OBJDIR:/=\/)$*.o/g' \
	| sed 's/$*.o/& $(subst /,\/,$@)/g' > $@

.PHONY : clean
clean:
	$(RMDIR) .obj/ .dep/
	$(RM) $(DYLIB)
	$(RM) $(FLENS_HOME)/$(DYLIB)
