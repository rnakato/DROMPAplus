/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_DRAW_H_
#define _DD_DRAW_H_

#include "dd_class.hpp"
#include "dd_readfile.hpp"
#include "WigStats.hpp"
#include "color.hpp"

class Figure {
  const chrsize &chr;
  vChrArray arrays;
  std::vector<SamplePairOverlayed> &vsamplepairoverlayed;
  std::vector<bed> regionBed;
  int32_t pagewidth;

  void setChIPInputRatio_SamplePairEach(const DROMPA::Global &p, SamplePairEach &pair, const std::string &chrname)
  {
    DEBUGprint("setSamplePairEachRatio");
    if (pair.argvInput == "") return;
    pair.ratio = 1;

    switch (p.getChIPInputNormType()) {
    case 0:  // not normalize
      pair.ratio = 1;
      break;
    case 1:  // total read for genome
      pair.ratio = getratio(arrays.arrays.at(pair.argvChIP).totalreadnum,
			    arrays.arrays.at(pair.argvInput).totalreadnum);
      break;
    case 2:  // total read for each chromosome
      pair.ratio = getratio(arrays.arrays.at(pair.argvChIP).totalreadnum_chr.at(chrname),
			    arrays.arrays.at(pair.argvInput).totalreadnum_chr.at(chrname));
      break;
    case 3:  // NCIS
      pair.ratio = 1;
      break;
    }
#ifdef DEBUG
    std::cout << "ChIP/Input Ratio for chr " << chrname << ": " << pair.ratio << std::endl;
#endif
  }

public:
  Figure(DROMPA::Global &p, const chrsize &_chr):
    chr(_chr),
    arrays(p, chr),
    vsamplepairoverlayed(p.samplepair),
    regionBed(p.drawregion.getRegionBedChr(chr.getname())),
    pagewidth(p.drawparam.width_draw_pixel)
  {
    for (auto &x: vsamplepairoverlayed) {
      setChIPInputRatio_SamplePairEach(p, x.first, chr.getname());
      if (x.OverlayExists()) setChIPInputRatio_SamplePairEach(p, x.second, chr.getname());
    }
  }

  int32_t Draw(DROMPA::Global &p) {
    if (p.drawregion.isRegionBed() && !regionBed.size()) return 0;
    std::cout << "Drawing.." << std::endl;
    DrawData(p);
    return 1;
  }

  void DrawData(DROMPA::Global &p) {
    DEBUGprint("Figure::DrawData");
    std::string pdffilename(p.getFigFileNameChr(chr.getrefname()));
    //  std::cout << chr.getrefname() << std::endl;

#ifdef CAIRO_HAS_PDF_SURFACE
    const auto surface = Cairo::PdfSurface::create(pdffilename, width, height);

    if (p.drawregion.isRegionBed()){  // --region
      int32_t region_no(1);
      for (auto &x: regionBed) {
	int32_t num_page(p.drawparam.getNumPage(x.start, x.end));
	for(int32_t i=0; i<num_page; ++i) {
	  std::cout << boost::format("   page %5d/%5d/%5d\r") % (i+1) % num_page % region_no << std::flush;
	  PDFPage page(p, arrays.arrays, vsamplepairoverlayed, surface, chr, x.start, x.end);
	  page.MakePage(p, i, std::to_string(region_no));
	}
	++region_no;
	printf("\n");
      }
    } else if (p.drawregion.isGeneLociFile()) {  // --genelocifile
      int32_t len(p.drawregion.getLenGeneLoci());
      for (auto &m: p.anno.gmp.at(rmchr(chr.getname()))) {
	if(!p.drawregion.ExistInGeneLociFile(m.second.gname)) continue;

	int32_t start = std::max(0, m.second.txStart - len);
	int32_t end   = std::min(m.second.txEnd + len, chr.getlen() -1);
	int32_t num_page(p.drawparam.getNumPage(start, end));
	for(int32_t i=0; i<num_page; ++i) {
	  std::cout << boost::format("   page %5d/%5d/%s\r") % (i+1) % num_page % m.second.gname << std::flush;
	  PDFPage page(p, arrays.arrays, vsamplepairoverlayed, surface, chr, start, end);
	  page.MakePage(p, i, m.second.gname);
	}
	printf("\n");
      }
    } else {  // whole chromosome
      int32_t num_page(p.drawparam.getNumPage(0, chr.getlen()));
      for (int32_t i=0; i<num_page; ++i) {
	std::cout << boost::format("   page %5d/%5d\r") % (i+1) % num_page << std::flush;
	PDFPage page(p, arrays.arrays, vsamplepairoverlayed, surface, chr, 0, chr.getlen());
	page.MakePage(p, i, "1");
      }
      printf("\n");
    }
    std::cout << "Wrote PDF file \"" << pdffilename << "\"" << std::endl;

#else
    std::cout << "You must compile cairo with PDF support for DROMPA+." << std::endl;
    return;
#endif
  }
};

#endif /* _DD_READFILE_H_ */
