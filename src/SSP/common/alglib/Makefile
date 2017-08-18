CC = g++
CFLAGS += -std=c++11 -O2 -Wall -W
LDFLAGS =
ALGLBDIR = alglib-3.10.0/src

.PHONY: all clean

all: libalglib.a

libalglib.a: alglib.cpp
	$(CC) -c $< $(ALGLBDIR)/specialfunctions.cpp $(ALGLBDIR)/ap.cpp $(ALGLBDIR)/alglibinternal.cpp $(CFLAGS)
	ar r $@ alglib.o alglibinternal.o ap.o specialfunctions.o
	rm alglib.o alglibinternal.o ap.o specialfunctions.o


clean:
	rm libalglib.a

libalglib.a: Makefile alglib.h
