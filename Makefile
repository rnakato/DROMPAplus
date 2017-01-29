CC = g++
CFLAGS += -std=c++11 -O2 -Wall -W
LDFLAGS =
LIBS += -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams
LIBS_DP += -lz -lgsl -lgslcblas -lboost_thread

SRCDIR = ./src
OBJDIR = ./obj
BINDIR = ./bin
SSPDIR = ./src/SSP
SSPSRCDIR = $(SSPDIR)/src
SSPCMNDIR = $(SSPDIR)/common
SSPOBJDIR = $(SSPDIR)/obj
SSPCMNOBJDIR = $(SSPDIR)/cobj
ALGLIBDIR = $(SSPDIR)/common/alglib

PROGRAMS = parse2wig+ drompa+
TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS))
#$(warning $(TARGET))

ifdef DEBUG
CFLAGS += -DDEBUG
endif

OBJS_UTIL = $(SSPCMNOBJDIR)/util.o $(SSPOBJDIR)/BoostOptions.o
OBJS_PW = $(OBJDIR)/pw_main.o $(SSPOBJDIR)/Mapfile.o $(SSPOBJDIR)/ParseMapfile.o $(OBJDIR)/pw_makefile.o $(OBJDIR)/ReadBpStatus.o $(SSPOBJDIR)/LibraryComplexity.o $(OBJDIR)/WigStats.o $(OBJDIR)/pw_gc.o $(SSPOBJDIR)/ssp_shiftprofile.o $(SSPCMNOBJDIR)/statistics.o $(ALGLIBDIR)/libalglib.a
OBJS_DD = $(OBJDIR)/dd_main.o $(OBJDIR)/dd_readfile.o $(OBJDIR)/ReadAnnotation.o

.PHONY: all clean

all: $(TARGET)

$(BINDIR)/parse2wig+: $(OBJS_PW) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS) $(LIBS_DP)
$(BINDIR)/drompa+: $(OBJS_DD) $(OBJS_UTIL)
	$(CC) -o $@ $^ $(LIBS) $(LIBS_DP)

$(ALGLIBDIR)/libalglib.a:
	make -C $(ALGLIBDIR)
$(SSPOBJDIR)/%.o: $(SSPSRCDIR)/%.cpp
	make -C $(SSPDIR) $(OBJDIR)/$(@F)
$(CMNOBJDIR)/%.o: $(CMNDIR)/%.cpp
	make -C $(SSPDIR) $(OBJDIR)/$(@F)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS) $(WFLAGS)

clean:
	rm -rf bin obj

HEADS_UTIL = $(SSPSRCDIR)/MThread.hpp $(SSPSRCDIR)/BoostOptions.hpp $(SSPCMNDIR)/inline.hpp $(SSPCMNDIR)/seq.hpp $(SSPCMNDIR)/util.hpp

$(OBJDIR)/dd_main.o: $(SRCDIR)/dd_opt.hpp
$(OBJDIR)/pw_main.o: $(SRCDIR)/pw_makefile.hpp $(SRCDIR)/pw_gc.hpp $(SRCDIR)/ReadBpStatus.hpp
$(OBJDIR)/pw_makefile.o: $(SRCDIR)/pw_makefile.hpp $(SRCDIR)/ReadBpStatus.hpp
$(OBJDIR)/pw_gc.o: $(SRCDIR)/pw_gc.hpp $(SRCDIR)/ReadBpStatus.hpp
$(OBJDIR)/ReadBpStatus.o: $(SRCDIR)/ReadBpStatus.hpp
$(OBJS_UTIL): Makefile $(HEADS_UTIL)
$(OBJS_PW): Makefile $(SRCDIR)/pw_gv.hpp $(SRCDIR)/WigStats.hpp $(SSPSRCDIR)/ParseMapfile.hpp $(SSPCMNDIR)/statistics.hpp $(SSPSRCDIR)/LibraryComplexity.hpp $(HEADS_UTIL)
$(OBJS_DD): Makefile $(SRCDIR)/dd_gv.hpp $(SRCDIR)/dd_readfile.hpp $(HEADS_UTIL)
