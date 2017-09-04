CC = g++
CFLAGS += -std=c++11 -O2 -Wall -W
LDFLAGS =
LIBS += -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams
LIBS_DP += -lz -lgsl -lgslcblas -lboost_thread
LIBS_CAIRO += `pkg-config gtkmm-3.0 --cflags --libs`

SRCDIR = ./src
OBJDIR = ./obj
BINDIR = ./bin
SSPDIR = ./src/SSP
SSPSRCDIR = $(SSPDIR)/src
SSPCMNDIR = $(SSPDIR)/common
SSPOBJDIR = $(SSPDIR)/obj
SSPCMNOBJDIR = $(SSPDIR)/cobj

PROGRAMS = parse2wig+ drompa+
TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS))
#$(warning $(TARGET))

ifdef DEBUG
CFLAGS += -DDEBUG
endif
ifdef PRINTFRAGMENT
CFLAGS += -DPRINTFRAGMENT
endif
ifdef PRINTREAD
CFLAGS += -DPRINTREAD
endif

OBJS_UTIL = $(SSPCMNOBJDIR)/util.o $(SSPCMNOBJDIR)/BoostOptions.o
OBJS_PW = $(OBJDIR)/pw_main.o $(OBJDIR)/pw_makefile.o $(OBJDIR)/WigStats.o $(OBJDIR)/GenomeCoverage.o $(OBJDIR)/GCnormalization.o 
OBJS_DD = $(OBJDIR)/dd_main.o $(OBJDIR)/dd_command.o $(OBJDIR)/dd_readfile.o $(OBJDIR)/dd_draw.o $(OBJDIR)/dd_peakcall.o $(OBJDIR)/WigStats.o  $(SSPCMNOBJDIR)/ReadAnnotation.o
OBJS_SSP = $(SSPOBJDIR)/Mapfile.o $(SSPOBJDIR)/ParseMapfile.o $(SSPOBJDIR)/ReadBpStatus.o $(SSPOBJDIR)/LibraryComplexity.o $(SSPOBJDIR)/ShiftProfile.o $(SSPCMNOBJDIR)/statistics.o $(SSPCMNOBJDIR)/ReadAnnotation.o $(SSPCMNOBJDIR)/util.o $(SSPCMNOBJDIR)/BoostOptions.o $(SSPCMNDIR)/BedFormat.hpp

.PHONY: all clean

all: $(TARGET)

$(BINDIR)/parse2wig+: $(OBJS_PW) $(OBJS_UTIL) $(OBJS_SSP)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS) $(LIBS_DP)  $(CFLAGS)

$(BINDIR)/drompa+: $(OBJS_DD) $(OBJS_UTIL) $(OBJS_SSP)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS) $(LIBS_DP) $(LIBS_CAIRO)  $(CFLAGS)

$(OBJDIR)/dd_draw.o: $(SRCDIR)/dd_draw.cpp
	$(CC) -o $@ -c $< $(CFLAGS) $(LIBS_CAIRO)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS) $(WFLAGS)

$(SSPOBJDIR)/%.o: $(SSPSRCDIR)/%.cpp
	$(MAKE) -C $(SSPDIR) $(OBJDIR)/$(notdir $@)

$(SSPCMNOBJDIR)/%.o: $(SSPCMNDIR)/%.cpp
	$(MAKE) -C $(SSPDIR) cobj/$(notdir $@)

clean:
	rm -rf bin obj
	make -C $(SSPDIR) clean

HEADS_UTIL = $(SSPSRCDIR)/MThread.hpp $(SSPCMNDIR)/BoostOptions.hpp $(SSPCMNDIR)/inline.hpp $(SSPCMNDIR)/seq.hpp $(SSPCMNDIR)/util.hpp

$(OBJDIR)/pw_main.o: $(SRCDIR)/pw_makefile.hpp $(SRCDIR)/GCnormalization.hpp $(SSPSRCDIR)/ReadBpStatus.hpp $(SRCDIR)/GenomeCoverage.hpp
$(OBJDIR)/pw_makefile.o: $(SRCDIR)/pw_makefile.hpp $(SSPSRCDIR)/ReadBpStatus.hpp
$(OBJDIR)/dd_main.o: $(SRCDIR)/dd_command.hpp
$(OBJDIR)/dd_command.o: $(SRCDIR)/dd_command.hpp
$(OBJDIR)/dd_draw.o: $(SRCDIR)/dd_draw_p.hpp
$(OBJDIR)/GCnormalization.o: $(SRCDIR)/GCnormalization.hpp $(SSPSRCDIR)/ReadBpStatus.hpp
$(OBJDIR)/ReadBpStatus.o: $(SSPSRCDIR)/ReadBpStatus.hpp
$(OBJS_UTIL): Makefile $(HEADS_UTIL)
$(OBJS_PW): Makefile $(SRCDIR)/pw_gv.hpp $(SRCDIR)/WigStats.hpp $(SSPSRCDIR)/ParseMapfile.hpp $(SSPCMNDIR)/statistics.hpp $(SSPSRCDIR)/LibraryComplexity.hpp $(HEADS_UTIL) $(SSPSRCDIR)/ShiftProfile.hpp
$(OBJS_DD): Makefile $(SRCDIR)/dd_readfile.hpp $(SRCDIR)/dd_class.hpp $(SRCDIR)/dd_draw.hpp $(SRCDIR)/dd_peakcall.hpp $(HEADS_UTIL)
