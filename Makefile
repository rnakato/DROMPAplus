CC = g++
CFLAGS += -std=c++11 -Wall -O2
LIBS += -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams
LIBS_DP += -lz -lgsl -lgslcblas -lboost_system -lboost_thread
TARGET = parse2wig+ drompa+ gtf2refFlat compare_bed2tss peak_occurance multibed2gene
ALGLBDIR = alglib-3.10.0/src

ifdef DEBUG
CFLAGS += -DDEBUG
endif

OBJS_UTIL = readdata.o util.o
OBJS_ALGLIB = $(ALGLBDIR)/specialfunctions.cpp $(ALGLBDIR)/ap.cpp $(ALGLBDIR)/alglibinternal.cpp
OBJS_GTF = gtf2refFlat.o
OBJS_COM = compare_bed2tss.o gene_bed.o
OBJS_PO = peak_occurance.o gene_bed.o alglib.o
OBJS_MG = multibed2gene.o gene_bed.o
OBJS_PW = pw_main.o pw_readmapfile.o pw_makefile.o pw_gc.o pw_shiftprofile.o statistics.o alglib.o
OBJS_DD = dd_main.o dd_readfile.o

HEADS_UTIL = util.h readdata.h macro.h seq.h
HEADS_PW = pw_gv.h pw_readmapfile.h pw_makefile.h pw_gc.h pw_shiftprofile.h statistics.h $(HEADS_UTIL)
HEADS_DD = dd_gv.h dd_readfile.h $(HEADS_UTIL)

.PHONY: all

all: $(TARGET)

echo:
	@echo "CFLAGS = $(CFLAGS)"

gtf2refFlat: $(OBJS_GTF) $(OBJS_UTIL)
	$(CC) -o $@ $^ $(LIBS)
compare_bed2tss: $(OBJS_COM) $(OBJS_UTIL)
	$(CC) -o $@ $^ $(LIBS)
peak_occurance: $(OBJS_PO) $(OBJS_UTIL)
	$(CC) -o $@ $^ $(OBJS_ALGLIB) $(LIBS)
multibed2gene: $(OBJS_MG) $(OBJS_UTIL)
	$(CC) -o $@ $^ $(LIBS)
parse2wig+: $(OBJS_PW) $(OBJS_UTIL)
	$(CC) -o $@ $^ $(OBJS_ALGLIB) $(LIBS) $(LIBS_DP)

drompa+: $(OBJS_DD) $(OBJS_UTIL)
	$(CC) -o $@ $^ $(LIBS) $(LIBS_DP)

.SUFFIXES: .o .cpp
.cpp.o:
	$(CC) -c $< $(CFLAGS)

clean:
	rm gtf2refFlat.o $(OBJS_COM) peak_occurance.o multibed2gene.o $(OBJS_PW) $(OBJS_DD) $(TARGET) *~

dd_main.o: dd_opt.h
pw_shiftprofile.o: pw_shiftprofile_p.h
$(OBJS_UTIL): Makefile $(HEADS_UTIL)
$(OBJS_PW): Makefile $(HEADS_PW) alglib.h
$(OBJS_DD): Makefile $(HEADS_DD)
$(OBJS_GTF): Makefile $(HEADS_UTIL)
$(OBJS_COM): Makefile $(HEADS_UTIL) gene_bed.h
$(OBJS_PO):  Makefile $(HEADS_UTIL) gene_bed.h alglib.h
$(OBJS_MG):  Makefile $(HEADS_UTIL) gene_bed.h
