CC = g++
CFLAGS += -std=c++11 -Wall -O2 -Iutil # -fopenmp
LIBS += -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lz -lgsl -lgslcblas -lboost_system -lboost_thread #-fopenmp 
TARGET = parse2wig+ drompa+
SRCDIR = util
ALGLBDIR = $(SRCDIR)/alglib-3.10.0/src

ifdef DEBUG
CFLAGS += -DDEBUG
endif

OBJS_UTIL = $(SRCDIR)/readdata.o $(SRCDIR)/warn.o
OBJS_ALGLIB = $(ALGLBDIR)/specialfunctions.cpp $(ALGLBDIR)/ap.cpp $(ALGLBDIR)/alglibinternal.cpp $(SRCDIR)/alglib.o
OBJS_COMMON = statistics.o util.o
OBJS_PW = pw_main.o pw_readmapfile.o pw_makefile.o pw_gc.o pw_shiftprofile.o $(OBJS_COMMON)
OBJS_DD = dd_main.o $(OBJS_COMMON)

HEADS_UTIL = common.h util.h statistics.h $(SRCDIR)/readdata.h $(SRCDIR)/warn.h $(SRCDIR)/macro.h $(SRCDIR)/seq.h
HEADS_PW = pw_gv.h pw_readmapfile.h pw_makefile.h pw_gc.h pw_shiftprofile.h pw_shiftprofile_p.h $(HEADS_UTIL)
HEADS_DD = dd_gv.h dd_opt.h $(HEADS_UTIL)

SUBDIRS := $(SRCDIR)
.PHONY: all $(SUBDIRS)

all: $(SUBDIRS) $(TARGET)

$(SUBDIRS):
	$(MAKE) -C $@

echo:
	@echo "CFLAGS = $(CFLAGS)"

parse2wig+: $(OBJS_PW)
	$(CC) -o $@ $(OBJS_UTIL) $^ $(OBJS_ALGLIB) $(LIBS)

drompa+: $(OBJS_DD)
	$(CC) -o $@ $(OBJS_UTIL) $^ $(LIBS)

pw_readmapfile.o: pw_readmapfile.cpp
	$(CC) -c $< $(CFLAGS)

.SUFFIXES: .o .cpp
.cpp.o:
	$(CC) -c $< $(CFLAGS)

clean:
	cd $(SRCDIR); make clean
	rm $(OBJS_PW) $(OBJS_DD) $(TARGET) *~

$(OBJS_PW): $(HEADS_PW) Makefile
$(OBJS_DD): $(HEADS_DD) Makefile
