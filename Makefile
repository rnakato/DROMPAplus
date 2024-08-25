.PHONY: all clean

BINDIR = bin
PROGRAMS = parse2wig+ drompa+
TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS))

ifdef DEBUG
CMAKEFLAGS = -DENABLE_DEBUG=ON
endif

SSPDIR = submodules/SSP
HTSLIBDIR = $(SSPDIR)/src/htslib-1.10.2/

all: $(TARGET) $(HTSLIBDIR)/libhts.a

$(TARGET): $(HTSLIBDIR)/libhts.a
	mkdir -p build && cd build && cmake $(CMAKEFLAGS) .. && make
	mkdir -p bin
	cp build/test/parse2wig/parse2wig+ build/test/drompa/drompa+ bin

$(HTSLIBDIR)/libhts.a:
	$(MAKE) -C $(HTSLIBDIR)

clean:
	rm -rf build bin
	make -C $(HTSLIBDIR) clean
	make -C $(SSPDIR) clean
