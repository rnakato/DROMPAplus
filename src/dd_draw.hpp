/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_DRAW_H_
#define _DD_DRAW_H_

#include "dd_gv.hpp"
#include "dd_readfile.hpp"

class Figure {
  vChrArray vReadArray;
  std::vector<SamplePairOverlayed> &vsamplepairoverlayed;
  std::vector<bed> regionBed;
  int32_t pagewidth;

public:
  Figure(DROMPA::Global &p, const chrsize &chr):
    vReadArray(p, chr),
    vsamplepairoverlayed(p.samplepair),
    regionBed(p.drawregion.getRegionBedChr(chr.getname())),
    pagewidth(p.drawparam.width_draw_pixel)
  {
    int32_t normtype(p.getChIPInputNormType());
    for (auto &x: vsamplepairoverlayed) {
      x.first.setScalingFactor(normtype, vReadArray, chr.getname());
      if (x.OverlayExists()) x.second.setScalingFactor(normtype, vReadArray, chr.getname());
    }
  }

  int32_t Draw(DROMPA::Global &p) {
    DEBUGprint_FUNCStart();

    if (p.drawregion.isRegionBed() && !regionBed.size()) return 0;
    std::cout << "Drawing.." << std::endl;

    std::string pdffilename(p.getFigFileNameChr(vReadArray.getchr().getrefname()));
    int32_t width(p.drawparam.width_page_pixel);
    int32_t height(p.drawparam.getPageHeight(p, vsamplepairoverlayed));

    if (p.drawregion.isRegionBed())         Draw_SpecificRegion(p, pdffilename, width, height);
    else if (p.drawregion.isGeneLociFile()) Draw_SpecificGene(p, pdffilename, width, height);
    else                                    Draw_WholeGenome(p, pdffilename, width, height);
    std::cout << "Wrote PDF file \"" << pdffilename << "\"" << std::endl;

    DEBUGprint_FUNCend();
    return 1;
  }

  void peakcall(DROMPA::Global &p, const std::string &chrname) {
    for (auto &x: vsamplepairoverlayed) {
      if (x.first.InputExists()) {
	x.first.peakcall_withInput(vReadArray, chrname, -log10(p.thre.pthre_inter), -log10(p.thre.pthre_enrich));
      }
      else x.first.peakcall_onlyChIP(vReadArray, chrname, -log10(p.thre.pthre_inter));
    }
  }

  void Draw_SpecificRegion(DROMPA::Global &p, std::string &pdffilename, int32_t width, int32_t height);
  void Draw_SpecificGene(DROMPA::Global &p, std::string &pdffilename, int32_t width, int32_t height);
  void Draw_WholeGenome(DROMPA::Global &p, std::string &pdffilename, int32_t width, int32_t height);
};

#endif /* _DD_DRAW_H_ */
