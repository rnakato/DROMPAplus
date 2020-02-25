CC = g++
CFLAGS = -std=c++11 -O2 -Wall -W -fuse-ld=gold
LDFLAGS =
LIBS += -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams
LIBS_DP += -lz -lgsl -lgslcblas -lboost_thread
LIBS_CAIRO += `pkg-config gtkmm-3.0 --cflags --libs`

SRCDIR = ./src
OBJDIR = ./obj
BINDIR = ./bin
SSPDIR = ./submodules/SSP
SSPSRCDIR = $(SSPDIR)/src
SSPCMNDIR = $(SSPDIR)/common
SSPOBJDIR = $(SSPDIR)/obj
SSPCMNOBJDIR = $(SSPDIR)/cobj

PROGRAMS = parse2wig+ drompa+
TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS))
#$(warning $(TARGET))

ifdef CLOCK
CFLAGS += -DCLOCK
OFLAGS += CLOCK=1
endif
ifdef DEBUG
CFLAGS += -DDEBUG
OFLAGS += DEBUG=1
endif
ifdef PRINTFRAGMENT
CFLAGS += -DPRINTFRAGMENT
endif
ifdef PRINTREAD
CFLAGS += -DPRINTREAD
endif

OBJS_UTIL = $(SSPCMNOBJDIR)/util.o $(SSPCMNOBJDIR)/BoostOptions.o
OBJS_PW = $(OBJDIR)/pw_main.o $(OBJDIR)/pw_makefile.o $(OBJDIR)/WigStats.o $(OBJDIR)/GenomeCoverage.o $(OBJDIR)/GCnormalization.o $(OBJDIR)/ReadAnnotation.o
OBJS_DD = $(OBJDIR)/dd_main.o $(OBJDIR)/dd_init.o $(OBJDIR)/dd_classfunc_draw.o $(OBJDIR)/dd_command.o $(OBJDIR)/dd_readfile.o $(OBJDIR)/dd_draw.o $(OBJDIR)/dd_chiadrop.o $(OBJDIR)/dd_drawgenes.o $(OBJDIR)/color.o $(OBJDIR)/dd_peakcall.o $(OBJDIR)/WigStats.o $(OBJDIR)/ReadAnnotation.o
OBJS_SSP = $(SSPOBJDIR)/Mapfile.o $(SSPOBJDIR)/ParseMapfile.o $(SSPOBJDIR)/ReadBpStatus.o $(SSPOBJDIR)/LibraryComplexity.o $(SSPOBJDIR)/ShiftProfile.o $(SSPCMNOBJDIR)/statistics.o $(SSPCMNOBJDIR)/util.o $(SSPCMNOBJDIR)/BoostOptions.o

.PHONY: all clean

all: $(TARGET) prnt

prnt: $(TARGET)
	@echo "\nAdd '$(CURDIR)/bin:$(CURDIR)/otherbins:$(CURDIR)/submodules/cpdf/Linux-Intel-64bit/' to your PATH."
#	@if [ ! -e scripts ]; then ln -s submodules/SSP/scripts/ .; fi

$(BINDIR)/parse2wig+: $(OBJS_PW) $(OBJS_UTIL) $(OBJS_SSP)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS) $(LIBS_DP) $(CFLAGS)

$(BINDIR)/drompa+: $(OBJS_DD) $(OBJS_UTIL) $(OBJS_SSP)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS) $(LIBS_DP) $(LIBS_CAIRO) $(CFLAGS)

$(OBJDIR)/dd_draw.o: $(SRCDIR)/dd_draw.cpp
	$(CC) -o $@ -c $< $(CFLAGS) $(LIBS_CAIRO)

$(OBJDIR)/dd_drawgenes.o: $(SRCDIR)/dd_drawgenes.cpp
	$(CC) -o $@ -c $< $(CFLAGS) $(LIBS_CAIRO)

$(OBJDIR)/dd_chiadrop.o: $(SRCDIR)/dd_chiadrop.cpp
	$(CC) -o $@ -c $< $(CFLAGS) $(LIBS_CAIRO)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS)

$(SSPOBJDIR)/%.o: $(SSPSRCDIR)/%.cpp
	$(MAKE) -C $(SSPDIR) $(OBJDIR)/$(notdir $@) $(OFLAGS)

$(SSPCMNOBJDIR)/%.o: $(SSPCMNDIR)/%.cpp
	$(MAKE) -C $(SSPDIR) cobj/$(notdir $@) $(OFLAGS)

clean:
	rm -rf bin obj
	make -C $(SSPDIR) clean

HEADS_UTIL = $(SSPSRCDIR)/MThread.hpp $(SSPCMNDIR)/BoostOptions.hpp $(SSPCMNDIR)/inline.hpp $(SSPCMNDIR)/seq.hpp $(SSPCMNDIR)/util.hpp $(SRCDIR)/WigStats.hpp $(SSPSRCDIR)/SeqStats.hpp

$(OBJDIR)/ReadAnnotation.o: $(SRCDIR)/ReadAnnotation.hpp $(SRCDIR)/GeneAnnotation.hpp
$(OBJDIR)/pw_main.o: $(SRCDIR)/pw_makefile.hpp $(SRCDIR)/GCnormalization.hpp $(SSPSRCDIR)/ReadBpStatus.hpp $(SRCDIR)/GenomeCoverage.hpp
$(OBJDIR)/pw_makefile.o: $(SRCDIR)/pw_makefile.hpp $(SSPSRCDIR)/ReadBpStatus.hpp
$(OBJDIR)/dd_main.o: $(SRCDIR)/dd_command.hpp $(SRCDIR)/dd_draw.hpp $(SRCDIR)/dd_profile.hpp
$(OBJDIR)/dd_command.o: $(SRCDIR)/dd_command.hpp
$(OBJDIR)/dd_classfunc_draw.o: $(SRCDIR)/dd_draw_environment_variable.hpp
$(OBJDIR)/dd_draw.o: $(SRCDIR)/dd_draw.hpp $(SRCDIR)/dd_draw_pdfpage.hpp $(SRCDIR)/dd_draw_myfunc.hpp $(SRCDIR)/dd_draw_environment_variable.hpp src/dd_draw_dataframe.hpp $(SRCDIR)/color.hpp
$(OBJDIR)/dd_drawgenes.o: $(SRCDIR)/dd_draw.hpp $(SRCDIR)/dd_draw_pdfpage.hpp $(SRCDIR)/dd_draw_myfunc.hpp $(SRCDIR)/dd_draw_environment_variable.hpp $(SRCDIR)/color.hpp
$(OBJDIR)/dd_chiadrop.o: $(SRCDIR)/dd_draw.hpp $(SRCDIR)/dd_draw_pdfpage.hpp $(SRCDIR)/dd_draw_myfunc.hpp $(SRCDIR)/dd_draw_environment_variable.hpp $(SRCDIR)/color.hpp
$(OBJDIR)/GCnormalization.o: $(SRCDIR)/GCnormalization.hpp $(SSPSRCDIR)/ReadBpStatus.hpp
$(OBJDIR)/ReadBpStatus.o: $(SSPSRCDIR)/ReadBpStatus.hpp
$(OBJS_UTIL): Makefile $(HEADS_UTIL)
$(OBJS_PW): Makefile $(SRCDIR)/version.hpp $(SRCDIR)/pw_gv.hpp $(SSPSRCDIR)/ParseMapfile.hpp $(SSPCMNDIR)/statistics.hpp $(SSPSRCDIR)/LibraryComplexity.hpp $(HEADS_UTIL) $(SSPSRCDIR)/ShiftProfile_p.hpp $(SSPSRCDIR)/ShiftProfile.hpp
$(OBJS_DD): Makefile $(SRCDIR)/version.hpp $(SRCDIR)/dd_readfile.hpp $(SRCDIR)/dd_gv.hpp $(SRCDIR)/dd_peakcall.hpp $(SSPCMNDIR)/BedFormat.hpp $(HEADS_UTIL) $(SRCDIR)/ReadAnnotation.hpp $(SRCDIR)/GeneAnnotation.hpp
