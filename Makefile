.PHONY: all build clean

BINDIR   = bin
BUILDDIR = build
PROGRAMS = parse2wig+ drompa+
TARGET   = $(addprefix $(BINDIR)/,$(PROGRAMS))

SSPDIR   = submodules/SSP
HTSLIBDIR = $(SSPDIR)/src/htslib-1.10.2

ifdef DEBUG
CMAKEFLAGS += -DENABLE_DEBUG=ON
endif

all: build

JOBS ?= 8

build: $(HTSLIBDIR)/libhts.a
	cmake -S . -B $(BUILDDIR) $(CMAKEFLAGS)
	cmake --build $(BUILDDIR) --parallel $(JOBS)
	mkdir -p $(BINDIR)
	cp $(BUILDDIR)/test/parse2wig/parse2wig+ $(BINDIR)/
	cp $(BUILDDIR)/test/drompa/drompa+ $(BINDIR)/


$(HTSLIBDIR)/libhts.a:
	$(MAKE) -C $(HTSLIBDIR)

clean:
	rm -rf $(BUILDDIR) $(BINDIR)
	$(MAKE) -C $(HTSLIBDIR) clean
	$(MAKE) -C $(SSPDIR) clean
