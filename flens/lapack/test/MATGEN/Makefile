TOPDIR = ..

include $(TOPDIR)/make.inc

SOURCE_FILES = $(wildcard *.f)
OBJECT_FILES = $(patsubst %.f, %.o, $(SOURCE_FILES))

all: $(TMGLIB)

$(TMGLIB) : $(OBJECT_FILES)
	ar cru $(TMGLIB) $(OBJECT_FILES)
	ranlib $(TMGLIB)


clean:
	$(RM) $(OBJECT_FILES) $(TMGLIB)
