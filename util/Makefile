CC = g++
CFLAGS += -std=c++11 -Wall -O2
TARGET = gtf2refFlat compare_bed2tss peak_occurance multibed2gene
SRCDIR = ./alglib-3.10.0/src

OBJS = warn.o readdata.o gene_bed.o
OBJSgtf = gtf2refFlat.o $(OBJS)
OBJScom = compare_bed2tss.o $(OBJS)
OBJSpo = peak_occurance.o alglib.o $(OBJS)
OBJSmg = multibed2gene.o $(OBJS)
OBALGLIB = $(SRCDIR)/specialfunctions.cpp $(SRCDIR)/ap.cpp $(SRCDIR)/alglibinternal.cpp

HEADS = warn.h readdata.h gene_bed.h cmdline.h seq.h macro.h

all: $(TARGET)

echo:
	@echo "CFLAGS = $(CFLAGS)"

gtf2refFlat: $(OBJSgtf)
	$(CC) -o $@ $^ 
compare_bed2tss: $(OBJScom)
	$(CC) -o $@ $^ 
peak_occurance: $(OBJSpo) 
	$(CC) -o $@ $^ $(OBALGLIB)
multibed2gene: $(OBJSmg)
	$(CC) -o $@ $^

.SUFFIXES: .o .cpp
.cpp.o:
	$(CC) -c $< $(CFLAGS)

clean:
	rm gtf2refFlat.o compare_bed2tss.o peak_occurance.o alglib.o multibed2gene.o $(OBJS) $(OBJS_alglib) $(TARGET)

$(OBJSgtf): $(HEADS) Makefile
$(OBJScom): $(HEADS) Makefile
$(OBJSpo):  $(HEADS) alglib.h Makefile
$(OBJSmg):  $(HEADS) Makefile
