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

  void DrawData(DROMPA::Global &p);
};

#endif /* _DD_READFILE_H_ */
