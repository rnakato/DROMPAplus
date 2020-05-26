.PHONY: all clean

SSPDIR = submodules/SSP
HTSLIBDIR = $(SSPDIR)/src/htslib-1.10.2/
SSPLIB = build/$(SSPDIR)/src/libssp_func.a build/$(SSPDIR)/common/libssp_common.a

all: bin/parse2wig+ bin/drompa+ $(HTSLIBDIR)/libhts.a

bin/parse2wig+: $(HTSLIBDIR)/libhts.a $(SSPLIB)
	mkdir -p build && cd build && cmake .. && make
	mkdir -p bin
	cp build/test/parse2wig/parse2wig+ bin

bin/drompa+: $(HTSLIBDIR)/libhts.a $(SSPLIB)
	mkdir -p build && cd build && cmake .. && make
	mkdir -p bin
	cp build/test/drompa/drompa+ bin

$(HTSLIBDIR)/libhts.a:
	$(MAKE) -C $(HTSLIBDIR)

clean:
	rm -rf build bin
	make -C $(HTSLIBDIR) clean
	make -C $(SSPDIR) clean
