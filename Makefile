CC = g++
CFLAGS += -std=c++11 -O2 -Wall -W
LDFLAGS =
LIBS += -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams
LIBS_DP += -lz -lgsl -lgslcblas -lboost_thread

SRCDIR = ./src
OBJDIR = ./obj
LIBDIR = ./lib
BINDIR = ./bin
SSPDIR = ./src/SSP
SSPSRCDIR = $(SSPDIR)/src
SSPOBJDIR = $(SSPDIR)/obj
ALGLIBDIR = $(SSPSRCDIR)/alglib

PROGRAMS = parse2wig+ drompa+
TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS))
#$(warning $(TARGET))

ifdef DEBUG
CFLAGS += -DDEBUG
endif

OBJS_UTIL = $(SSPOBJDIR)/readdata.o $(SSPOBJDIR)/util.o
OBJS_PW = $(OBJDIR)/pw_main.o $(SSPOBJDIR)/Mapfile.o $(SSPOBJDIR)/pw_readmapfile.o $(OBJDIR)/pw_makefile.o $(OBJDIR)/readbpstatus.o $(SSPOBJDIR)/LibraryComplexity.o $(OBJDIR)/pw_gc.o $(SSPOBJDIR)/ssp_shiftprofile.o $(SSPOBJDIR)/statistics.o $(ALGLIBDIR)/libalglib.a
OBJS_DD = $(OBJDIR)/dd_main.o $(OBJDIR)/dd_readfile.o

.PHONY: all clean

all: $(TARGET)

$(BINDIR)/parse2wig+: $(OBJS_PW) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS) $(LIBS_DP)
$(BINDIR)/drompa+: $(OBJS_DD) $(OBJS_UTIL)
	$(CC) -o $@ $^ $(LIBS) $(LIBS_DP)

$(LIBDIR)/libalglib.a:
	make -C $(ALGLIBDIR)
$(SSPOBJDIR)/%.o: $(SSPSRCDIR)/%.cpp
	make -C $(SSPDIR) $(OBJDIR)/$(@F)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS) $(WFLAGS)

clean:
	rm -rf bin obj

HEADS_UTIL = $(SSPSRCDIR)/util.h $(SSPSRCDIR)/readdata.h $(SSPSRCDIR)/macro.h $(SSPSRCDIR)/seq.h $(SSPSRCDIR)/mthread.h

$(OBJDIR)/dd_main.o: $(SRCDIR)/dd_opt.h
$(OBJDIR)/pw_main.o: $(SRCDIR)/pw_makefile.h $(SRCDIR)/pw_gc.h $(SRCDIR)/readbpstatus.h
$(OBJDIR)/pw_readmapfile.o: $(SSPSRCDIR)/ssp_shiftprofile.h
$(OBJDIR)/pw_makefile.o: $(SRCDIR)/pw_makefile.h $(SRCDIR)/readbpstatus.h
$(OBJDIR)/pw_gc.o: $(SRCDIR)/pw_gc.h $(SRCDIR)/readbpstatus.h
$(OBJDIR)/readbpstatus.o: $(SRCDIR)/readbpstatus.h
$(OBJS_UTIL): Makefile $(HEADS_UTIL)
$(OBJS_PW): Makefile $(SSPSRCDIR)/pw_gv.h $(SSPSRCDIR)/pw_readmapfile.h $(SSPSRCDIR)/statistics.h $(SSPSRCDIR)/LibraryComplexity.hpp $(HEADS_UTIL)
$(OBJS_DD): Makefile $(SRCDIR)/dd_gv.h $(SRCDIR)/dd_readfile.h $(HEADS_UTIL)
