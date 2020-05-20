CC = clang++
CFLAGS = -std=c++11 -O2 -Wall -W
LIBS += -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lpthread
LIBS_DP += -lz -lgsl -lgslcblas -lboost_thread -lcurl -llzma -lbz2
LIBS_CAIRO = `pkg-config gtkmm-3.0 --cflags --silence-errors` `pkg-config gtkmm-2.4 --cflags --silence-errors`
FLAG_CAIRO = `pkg-config gtkmm-3.0 --libs --silence-errors` `pkg-config gtkmm-2.4 --libs --silence-errors`

SRCDIR = ./src
OBJDIR = ./obj
BINDIR = ./bin
SSPDIR = ./submodules/SSP
SSPSRCDIR = $(SSPDIR)/src
SSPCMNDIR = $(SSPDIR)/common
SSPOBJDIR = $(SSPDIR)/obj
SSPCMNOBJDIR = $(SSPDIR)/cobj
HTSLIBDIR = $(SSPDIR)/src/htslib-1.10.2

PROGRAMS = parse2wig+ drompa+
TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS))

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

OBJS_UTIL = $(SSPCMNOBJDIR)/util.o $(SSPCMNOBJDIR)/BoostOptions.o $(SSPCMNOBJDIR)/statistics.o $(SSPCMNOBJDIR)/gzstream.o
OBJS_PW = $(OBJDIR)/pw_main.o $(OBJDIR)/pw_makefile.o $(OBJDIR)/WigStats.o $(OBJDIR)/GenomeCoverage.o $(OBJDIR)/GCnormalization.o $(OBJDIR)/ReadAnnotation.o $(OBJDIR)/ReadMpbldata.o $(OBJDIR)/significancetest.o $(OBJDIR)/pw_strShiftProfile.o
OBJS_DD = $(OBJDIR)/dd_main.o $(OBJDIR)/dd_init.o $(OBJDIR)/dd_draw_dataframe.o $(OBJDIR)/dd_classfunc_draw.o $(OBJDIR)/dd_command.o $(OBJDIR)/dd_readfile.o $(OBJDIR)/dd_draw.o $(OBJDIR)/dd_chiadrop.o $(OBJDIR)/dd_drawgenes.o $(OBJDIR)/color.o $(OBJDIR)/significancetest.o $(OBJDIR)/WigStats.o $(OBJDIR)/ReadAnnotation.o $(OBJDIR)/dd_sample_definition.o $(OBJDIR)/dd_profile.o
OBJS_SSP = $(SSPOBJDIR)/Mapfile.o $(SSPOBJDIR)/ParseMapfile.o $(SSPOBJDIR)/LibraryComplexity.o $(SSPOBJDIR)/ShiftProfile.o

.PHONY: all clean

all: $(TARGET) prnt

prnt: $(TARGET)
	@echo "\nAdd '$(CURDIR)/bin:$(CURDIR)/otherbins:$(CURDIR)/submodules/cpdf/Linux-Intel-64bit/' to your PATH."

$(BINDIR)/parse2wig+: $(OBJS_PW) $(OBJS_UTIL) $(OBJS_SSP) $(HTSLIBDIR)/libhts.a
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS) $(LIBS_DP) $(CFLAGS)

$(BINDIR)/drompa+: $(OBJS_DD) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS) $(LIBS_DP) $(FLAG_CAIRO) $(CFLAGS)

$(OBJDIR)/dd_draw.o: $(SRCDIR)/dd_draw.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS) $(LIBS_CAIRO)

$(OBJDIR)/dd_draw_dataframe.o: $(SRCDIR)/dd_draw_dataframe.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS) $(LIBS_CAIRO)

$(OBJDIR)/dd_drawgenes.o: $(SRCDIR)/dd_drawgenes.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS) $(LIBS_CAIRO)

$(OBJDIR)/dd_chiadrop.o: $(SRCDIR)/dd_chiadrop.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS) $(LIBS_CAIRO)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS)

$(SSPOBJDIR)/%.o: $(SSPSRCDIR)/%.cpp
	$(MAKE) -C $(SSPDIR) $(OBJDIR)/$(notdir $@) $(OFLAGS)

$(SSPCMNOBJDIR)/gzstream.o: $(SSPCMNDIR)/gzstream.C $(SSPCMNDIR)/gzstream.h
	$(MAKE) -C $(SSPDIR) cobj/gzstream.o

$(SSPCMNOBJDIR)/%.o: $(SSPCMNDIR)/%.cpp
	$(MAKE) -C $(SSPDIR) cobj/$(notdir $@) $(OFLAGS)

$(HTSLIBDIR)/libhts.a:
	$(MAKE) -C $(HTSLIBDIR)

clean:
	rm -rf bin obj
	make -C $(SSPDIR) clean

HEADS_UTIL = $(SSPSRCDIR)/MThread.hpp $(SSPCMNDIR)/BoostOptions.hpp $(SSPCMNDIR)/inline.hpp $(SSPCMNDIR)/seq.hpp $(SSPCMNDIR)/util.hpp $(SRCDIR)/WigStats.hpp $(SSPSRCDIR)/SeqStats.hpp

$(OBJDIR)/ReadAnnotation.o: $(SRCDIR)/ReadAnnotation.hpp $(SRCDIR)/GeneAnnotation.hpp
$(OBJDIR)/ReadMpbldata.o: $(SRCDIR)/ReadMpbldata.hpp $(SRCDIR)/BpStatus.hpp
$(OBJDIR)/pw_main.o: $(SRCDIR)/pw_makefile.hpp $(SRCDIR)/GCnormalization.hpp $(SRCDIR)/GenomeCoverage.hpp
$(OBJDIR)/pw_makefile.o: $(SRCDIR)/pw_makefile.hpp $(SRCDIR)/ReadMpbldata.hpp
$(OBJDIR)/dd_main.o: $(SRCDIR)/dd_command.hpp $(SRCDIR)/dd_draw.hpp $(SRCDIR)/dd_profile.hpp
$(OBJDIR)/dd_command.o: $(SRCDIR)/dd_command.hpp
$(OBJDIR)/dd_classfunc_draw.o: $(SRCDIR)/dd_draw_environment_variable.hpp
$(OBJDIR)/dd_draw.o: $(SRCDIR)/dd_draw.hpp $(SRCDIR)/dd_draw_pdfpage.hpp $(SRCDIR)/dd_draw_myfunc.hpp $(SRCDIR)/dd_draw_environment_variable.hpp $(SRCDIR)/dd_draw_dataframe.hpp $(SRCDIR)/color.hpp
$(OBJDIR)/dd_drawgenes.o: $(SRCDIR)/dd_draw.hpp $(SRCDIR)/dd_draw_pdfpage.hpp $(SRCDIR)/dd_draw_myfunc.hpp $(SRCDIR)/dd_draw_environment_variable.hpp $(SRCDIR)/color.hpp $(SRCDIR)/dd_draw_dataframe.hpp
$(OBJDIR)/dd_chiadrop.o: $(SRCDIR)/dd_draw.hpp $(SRCDIR)/dd_draw_pdfpage.hpp $(SRCDIR)/dd_draw_myfunc.hpp $(SRCDIR)/dd_draw_environment_variable.hpp $(SRCDIR)/color.hpp $(SRCDIR)/dd_draw_dataframe.hpp
$(OBJDIR)/dd_draw_dataframe.o: $(SRCDIR)/dd_draw.hpp $(SRCDIR)/dd_draw_pdfpage.hpp $(SRCDIR)/dd_draw_myfunc.hpp $(SRCDIR)/dd_draw_environment_variable.hpp $(SRCDIR)/color.hpp $(SRCDIR)/dd_draw_dataframe.hpp
$(OBJDIR)/GCnormalization.o: $(SRCDIR)/GCnormalization.hpp $(SRCDIR)/ReadMpbldata.hpp
$(OBJS_UTIL): Makefile $(HEADS_UTIL)
$(OBJS_PW): Makefile $(SRCDIR)/version.hpp $(SRCDIR)/pw_gv.hpp $(SSPCMNDIR)/statistics.hpp $(SSPSRCDIR)/LibraryComplexity.hpp $(HEADS_UTIL) $(SSPSRCDIR)/ShiftProfile_p.hpp $(SSPSRCDIR)/ShiftProfile.hpp $(SRCDIR)/significancetest.hpp $(SRCDIR)/BpStatus.hpp $(SRCDIR)/pw_strShiftProfile.hpp $(SRCDIR)/SeqStatsDROMPA.hpp
$(OBJS_DD): Makefile $(SRCDIR)/version.hpp $(SRCDIR)/dd_readfile.hpp $(SRCDIR)/dd_gv.hpp $(SRCDIR)/extendBedFormat.hpp $(SRCDIR)/dd_sample_definition.hpp $(SRCDIR)/significancetest.hpp $(HEADS_UTIL) $(SRCDIR)/ReadAnnotation.hpp $(SRCDIR)/GeneAnnotation.hpp
