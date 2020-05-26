.PHONY: all clean

HTSLIBDIR=submodules/SSP/src/htslib-1.10.2/

all: bin/parse2wig+ bin/drompa+ $(HTSLIBDIR)/libhts.a

bin/parse2wig+: $(HTSLIBDIR)/libhts.a
	mkdir -p build && cd build && cmake .. && make
	mkdir -p bin
	cp build/test/parse2wig+ build/test/drompa+ bin

bin/drompa+: $(HTSLIBDIR)/libhts.a
	mkdir -p build && cd build && cmake .. && make
	mkdir -p bin
	cp build/test/parse2wig/parse2wig+ build/test/drompa/drompa+ bin

$(HTSLIBDIR)/libhts.a:
	$(MAKE) -C $(HTSLIBDIR)

clean:
	rm -rf build bin
	make -C $(HTSLIBDIR) clean
