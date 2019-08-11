/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_DRAW_H_
#define _DD_DRAW_H_

#include "dd_gv.hpp"
#include "dd_readfile.hpp"
#include "WigStats.hpp"
#include "color.hpp"

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
      x.first.setratio(normtype, vReadArray, chr.getname());
      if (x.OverlayExists()) x.second.setratio(normtype, vReadArray, chr.getname());
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
    std::string pdffilename(p.getFigFileNameChr(vReadArray.getchr().getrefname()));
    int32_t width(p.drawparam.width_page_pixel);
    int32_t height(p.drawparam.getPageHeight(p, vsamplepairoverlayed));

    if (p.drawregion.isRegionBed()) DrawData_Region(p, pdffilename, width, height);
    else if (p.drawregion.isGeneLociFile()) DrawData_GeneLoci(p, pdffilename, width, height);
    else  DrawData_Whole(p, pdffilename, width, height);

    std::cout << "Wrote PDF file \"" << pdffilename << "\"" << std::endl;
  }

  void DrawData_Region(DROMPA::Global &p, std::string &pdffilename, int32_t width, int32_t height);
  void DrawData_GeneLoci(DROMPA::Global &p, std::string &pdffilename, int32_t width, int32_t height);
  void DrawData_Whole(DROMPA::Global &p, std::string &pdffilename, int32_t width, int32_t height);
};

#endif /* _DD_DRAW_H_ */
