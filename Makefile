CC = g++
CFLAGS += -std=c++11 -Wall -O2 -Iutil 
LIBS += -lboost_program_options -lboost_system -lboost_filesystem -lz
TARGET = parse2wig+ drompa+
SRCDIR = util

ifdef DEBUG
CFLAGS += -DDEBUG
endif

OBJS_UTIL = $(SRCDIR)/readdata.o $(SRCDIR)/warn.o
OBJS_COMMON = util.o
OBJS_PW = pw_main.o pw_readmapfile.o pw_makefile.o pw_gc.o $(OBJS_COMMON)
OBJS_DD = drompa+.o $(OBJS_COMMON)

HEADS_UTIL = common.h util.h $(SRCDIR)/readdata.h $(SRCDIR)/warn.h $(SRCDIR)/macro.h $(SRCDIR)/seq.h
HEADS_PW = pw_gv.h pw_readmapfile.h pw_makefile.h pw_gc.h $(HEADS_UTIL)
HEADS_DD = dd_gv.h dd_opt.h $(HEADS_UTIL)

SUBDIRS := $(SRCDIR)
.PHONY: all $(SUBDIRS)

all: $(SUBDIRS) $(TARGET)

$(SUBDIRS):
	$(MAKE) -C $@

echo:
	@echo "CFLAGS = $(CFLAGS)"

parse2wig+: $(OBJS_PW)
	$(CC) -o $@ $(OBJS_UTIL) $^ $(LIBS)

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
