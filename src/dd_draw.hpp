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

  void DrawData(DROMPA::Global &p);
};

#endif /* _DD_DRAW_H_ */
