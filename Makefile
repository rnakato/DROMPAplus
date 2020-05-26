.PHONY: all clean

BINDIR = bin
PROGRAMS = parse2wig+ drompa+
TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS))

SSPDIR = submodules/SSP
HTSLIBDIR = $(SSPDIR)/src/htslib-1.10.2/
#SSPLIB = build/$(SSPDIR)/src/libssp_func.a build/$(SSPDIR)/common/libssp_common.a

all: $(TARGET) $(HTSLIBDIR)/libhts.a

$(TARGET): $(HTSLIBDIR)/libhts.a
	mkdir -p build && cd build && cmake .. && make
	mkdir -p bin
	cp build/test/parse2wig/parse2wig+ build/test/drompa/drompa+ bin

$(HTSLIBDIR)/libhts.a:
	$(MAKE) -C $(HTSLIBDIR)

clean:
	rm -rf build bin
	make -C $(HTSLIBDIR) clean
	make -C $(SSPDIR) clean
