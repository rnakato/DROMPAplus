CC = g++
CFLAGS += -std=c++11 -Wall -O2
LIBS += -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lz -lgsl -lgslcblas -lboost_system -lboost_thread
TARGET = parse2wig+ drompa+ gtf2refFlat compare_bed2tss peak_occurance multibed2gene
ALGLBDIR = alglib-3.10.0/src

ifdef DEBUG
CFLAGS += -DDEBUG
endif

OBJS_UTIL = readdata.o statistics.o util.o
OBJS_ALGLIB = $(ALGLBDIR)/specialfunctions.cpp $(ALGLBDIR)/ap.cpp $(ALGLBDIR)/alglibinternal.cpp
OBJS_GTF = gtf2refFlat.o $(OBJS_UTIL)
OBJS_COM = compare_bed2tss.o gene_bed.o $(OBJS_UTIL)
OBJS_PO = peak_occurance.o gene_bed.o alglib.o $(OBJS_UTIL)
OBJS_MG = multibed2gene.o gene_bed.o $(OBJS_UTIL)
OBJS_PW = pw_main.o pw_readmapfile.o pw_makefile.o pw_gc.o pw_shiftprofile.o alglib.o $(OBJS_UTIL)
OBJS_DD = dd_main.o dd_readfile.o $(OBJS_UTIL)

HEADS_UTIL = common.h util.h statistics.h readdata.h macro.h seq.h gene_bed.h
HEADS_PW = pw_gv.h pw_readmapfile.h pw_makefile.h pw_gc.h pw_shiftprofile.h $(HEADS_UTIL)
HEADS_DD = dd_gv.h dd_readfile.h $(HEADS_UTIL)

.PHONY: all

all: $(TARGET)

echo:
	@echo "CFLAGS = $(CFLAGS)"

gtf2refFlat: $(OBJS_GTF)
	$(CC) -o $@ $^ $(LIBS)
compare_bed2tss: $(OBJS_COM)
	$(CC) -o $@ $^ $(LIBS)
peak_occurance: $(OBJS_PO)
	$(CC) -o $@ $^ $(OBJS_ALGLIB) $(LIBS)
multibed2gene: $(OBJS_MG)
	$(CC) -o $@ $^ $(LIBS)
parse2wig+: $(OBJS_PW)
	$(CC) -o $@ $^ $(OBJS_ALGLIB) $(LIBS)

drompa+: $(OBJS_DD)
	$(CC) -o $@ $^ $(LIBS)

.SUFFIXES: .o .cpp
.cpp.o:
	$(CC) -c $< $(CFLAGS)

clean:
	rm gtf2refFlat.o $(OBJS_COM) $(OBJS_PO) multibed2gene.o $(OBJS_PW) $(OBJS_DD) $(OBJS_COMMON) $(TARGET) *~

dd_main.o: dd_opt.h
pw_shiftprofile.o: pw_shiftprofile_p.h
$(OBJS_COMMON): $(HEADS_UTIL)
$(OBJS_PW): $(HEADS_PW) Makefile
$(OBJS_DD): $(HEADS_DD) Makefile
$(OBJSgtf): $(HEADS) Makefile
$(OBJScom): $(HEADS) Makefile
$(OBJSpo):  $(HEADS) alglib.h Makefile
$(OBJSmg):  $(HEADS) Makefile
