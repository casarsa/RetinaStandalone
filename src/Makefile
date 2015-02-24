ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -fPIC -ansi -D_GNU_SOURCE -O2 -Wall -Wextra 
LDFLAGS       = -g3 

SOFLAGS       = -shared

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)

CXXFLAGS     += $(ROOTCFLAGS)
NGLIBS        = $(ROOTGLIBS) 
GLIBS         = $(filter-out -lNew, $(NGLIBS))

OBJDIR        = ../tmp
INCLUDEDIR    = ../interface
CXX	     += -I$(INCLUDEDIR)


HEADERS = Retina.h Track.h 
DEPS = $(patsubst %,$(INCLUDEDIR)/%,$(HEADERS))

OBJFILES = RetinaStandalone.o Retina.o Track.o 
OBJ = $(patsubst %,$(OBJDIR)/%,$(OBJFILES))

###########
# TARGETS #
###########

all: Retina

$(OBJDIR)/%.o: %.cc $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

Retina: $(OBJ)
	$(CXX) -o ../bin/$@ $^  $(LDFLAGS) $(GLIBS)	

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o *~ core $(INCLUDEDIR)/*~ 


