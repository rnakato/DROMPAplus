CC = g++
CFLAGS += -std=c++11 -O2 -Wall -W
LDFLAGS =
LIBS += -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams
LIBS_DP += -lz -lgsl -lgslcblas -lboost_thread

SRCDIR = ./src
OBJDIR = ./obj
LIBDIR = ./lib
BINDIR = ./bin
ALGLBDIR = $(SRCDIR)/alglib-3.10.0/src

PROGRAMS = parse2wig+ drompa+ gtf2refFlat compare_bed2tss peak_occurance multibed2gene
TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS))
#$(warning $(TARGET))

ifdef DEBUG
CFLAGS += -DDEBUG
endif

OBJS_UTIL = $(OBJDIR)/readdata.o $(OBJDIR)/util.o
OBJS_GTF = $(OBJDIR)/gtf2refFlat.o
OBJS_COM = $(OBJDIR)/compare_bed2tss.o $(OBJDIR)/gene_bed.o
OBJS_PO = $(OBJDIR)/peak_occurance.o $(OBJDIR)/gene_bed.o $(LIBDIR)/libalglib.a
OBJS_MG = $(OBJDIR)/multibed2gene.o $(OBJDIR)/gene_bed.o
OBJS_PW = $(OBJDIR)/pw_main.o $(OBJDIR)/pw_readmapfile.o $(OBJDIR)/pw_makefile.o $(OBJDIR)/pw_gc.o $(OBJDIR)/pw_shiftprofile.o $(OBJDIR)/statistics.o $(LIBDIR)/libalglib.a
OBJS_DD = $(OBJDIR)/dd_main.o $(OBJDIR)/dd_readfile.o

.PHONY: all clean

all: $(TARGET)

$(BINDIR)/gtf2refFlat: $(OBJS_GTF) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS)
$(BINDIR)/compare_bed2tss: $(OBJS_COM) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS)
$(BINDIR)/peak_occurance: $(OBJS_PO) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS)
$(BINDIR)/multibed2gene: $(OBJS_MG) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS)
$(BINDIR)/parse2wig+: $(OBJS_PW) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS) $(LIBS_DP)
$(BINDIR)/drompa+: $(OBJS_DD) $(OBJS_UTIL)
	$(CC) -o $@ $^ $(LIBS) $(LIBS_DP)

$(LIBDIR)/libalglib.a: $(SRCDIR)/alglib.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -c $< $(ALGLBDIR)/specialfunctions.cpp $(ALGLBDIR)/ap.cpp $(ALGLBDIR)/alglibinternal.cpp $(CFLAGS)
	ar r $@ alglib.o alglibinternal.o ap.o specialfunctions.o
	rm alglib.o alglibinternal.o ap.o specialfunctions.o

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS) $(WFLAGS)

clean:
	rm -rf bin lib obj

HEADS_UTIL = $(SRCDIR)/util.h $(SRCDIR)/readdata.h $(SRCDIR)/macro.h $(SRCDIR)/seq.h

$(SRCDIR)/alglib.o: Makefile $(SRCDIR)/alglib.h
$(OBJDIR)/dd_main.o: $(SRCDIR)/dd_opt.h
$(OBJDIR)/pw_main.o: $(SRCDIR)/pw_makefile.h $(SRCDIR)/pw_gc.h
$(OBJDIR)/pw_readmapfile.o: $(SRCDIR)/pw_shiftprofile.h
$(OBJDIR)/pw_makefile.o: $(SRCDIR)/pw_makefile.h
$(OBJDIR)/pw_gc.o: $(SRCDIR)/pw_gc.h
$(OBJDIR)/pw_shiftprofile.o: Makefile $(SRCDIR)/pw_shiftprofile_p.h $(SRCDIR)/pw_shiftprofile.h
$(OBJS_UTIL): Makefile $(HEADS_UTIL)
$(OBJS_PW): Makefile $(SRCDIR)/pw_gv.h $(SRCDIR)/pw_readmapfile.h $(SRCDIR)/statistics.h $(HEADS_UTIL)
$(OBJS_DD): Makefile $(SRCDIR)/dd_gv.h $(SRCDIR)/dd_readfile.h $(HEADS_UTIL)
$(OBJS_GTF): Makefile $(HEADS_UTIL)
$(OBJS_COM): Makefile $(HEADS_UTIL) $(SRCDIR)/gene_bed.h
$(OBJS_PO):  Makefile $(HEADS_UTIL) $(SRCDIR)/gene_bed.h
$(OBJS_MG):  Makefile $(HEADS_UTIL) $(SRCDIR)/gene_bed.h
