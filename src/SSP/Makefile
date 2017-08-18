CC = g++
CFLAGS += -std=c++11 -O2 -Wall -W
LDFLAGS = -lz -lgsl -lgslcblas -lboost_thread
LIBS += -lboost_program_options -lboost_system -lboost_filesystem 

SRCDIR = ./src
CMNDIR = ./common
OBJDIR = ./obj
CMNOBJDIR = ./cobj
ALGLIBDIR = ./common/alglib
BINDIR = ./bin

PROGRAMS = ssp 
TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS))
#$(warning $(TARGET))

ifdef DEBUG
CFLAGS += -DDEBUG
endif
ifdef PRINTREAD
CFLAGS += -DPRINTREAD
endif

OBJS = $(OBJDIR)/ssp_main.o $(OBJDIR)/Mapfile.o $(OBJDIR)/ParseMapfile.o $(OBJDIR)/ReadBpStatus.o $(OBJDIR)/LibraryComplexity.o $(OBJDIR)/ShiftProfile.o $(OBJDIR)/FragmentClusterScore.o
OBJS += $(CMNOBJDIR)/statistics.o $(CMNOBJDIR)/ReadAnnotation.o $(CMNOBJDIR)/util.o $(CMNOBJDIR)/BoostOptions.o $(ALGLIBDIR)/libalglib.a

.PHONY: all clean

all: $(TARGET)

$(BINDIR)/ssp: $(OBJS)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS) $(LDFLAGS)

$(ALGLIBDIR)/libalglib.a:
	make -C $(ALGLIBDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS) $(WFLAGS)

$(CMNOBJDIR)/%.o: $(CMNDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS) $(WFLAGS)

clean:
	rm -rf $(BINDIR) $(OBJDIR) $(CMNOBJDIR)
	make clean -C $(ALGLIBDIR)

HEADS = $(SRCDIR)/ssp_gv.hpp $(SRCDIR)/Mapfile.hpp $(SRCDIR)/ParseMapfile.hpp $(SRCDIR)/LibraryComplexity.hpp $(CMNDIR)/BoostOptions.hpp $(SRCDIR)/MThread.hpp $(SRCDIR)/SeqStats.hpp $(SRCDIR)/BpStatus.hpp $(CMNDIR)/BedFormat.hpp $(SRCDIR)/ReadBpStatus.hpp $(SRCDIR)/FragmentClusterScore.hpp
HEADS += $(CMNDIR)/inline.hpp $(CMNDIR)/seq.hpp $(CMNDIR)/statistics.hpp $(CMNDIR)/util.hpp

$(OBJDIR)/ParseMapfile.o: $(SRCDIR)/ShiftProfile.hpp
$(OBJDIR)/ShiftProfile.o: $(SRCDIR)/ShiftProfile_p.hpp $(SRCDIR)/ShiftProfile.hpp
$(OBJDIR)/FragmentCluterScore.o: $(SRCDIR)/FragmentCluterScore_p.hpp
$(OBJS): Makefile $(HEADS)
