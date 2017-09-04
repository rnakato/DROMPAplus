/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_DRAW_H_
#define _DD_DRAW_H_

#include "dd_class.hpp"
#include "dd_readfile.hpp"
#include "WigStats.hpp"

class ChrArray {
public:
  int32_t binsize;
  int32_t nbin;
  WigArray array;
  WigStats stats;
  int32_t totalreadnum;
  std::unordered_map<std::string, int32_t> totalreadnum_chr;

  ChrArray(){}
  ChrArray(const DROMPA::Global &p,
	   const std::pair<const std::string, SampleFile> &x,
	   const chrsize &chr):
    binsize(x.second.getbinsize()), nbin(chr.getlen()/binsize +1),
    array(loadWigData(x.first, x.second, chr)),
    stats(nbin, p.thre.pthre_inter),
    totalreadnum(x.second.gettotalreadnum()),
    totalreadnum_chr(x.second.gettotalreadnum_chr())
  {
    if(p.getSmoothing()) array.Smoothing(p.getSmoothing());
    stats.setWigStats(array);
    stats.peakcall(array, chr.getname());
  }
};
  
class SamplePairChr {
  int32_t binsize;
public:
  std::string argvChIP;
  std::string argvInput;
  std::string label;
  double ratio;
  bool isPeaks;
  std::vector<bed> peaks;

  SamplePairChr(const DROMPA::Global &p, const SamplePair &x,
		const std::unordered_map<std::string, ChrArray> &arrays,
		const std::string &chrname):
    binsize(arrays.at(x.argvChIP).binsize),
    argvChIP(x.argvChIP), argvInput(x.argvInput),
    label(x.label), ratio(0), isPeaks(false)
  {
    if (argvInput != "") setRatio(p, arrays, chrname);
    if (x.peak_argv != "") {
      isPeaks = true;
      peaks = x.getpeaksChr(chrname);
    }
  }
  
  void setRatio (const DROMPA::Global &p, const std::unordered_map<std::string, ChrArray> &arrays, const std::string &chrname) {
    DEBUGprint("setRatio");
    switch(p.getNorm()) {
    case 0:
      ratio = 1;
      break;
    case 1:
      ratio = getratio(arrays.at(argvChIP).totalreadnum, arrays.at(argvInput).totalreadnum);
      break;
    case 2:
      ratio = getratio(arrays.at(argvChIP).totalreadnum_chr.at(chrname), arrays.at(argvInput).totalreadnum_chr.at(chrname));
      break;
    case 3:
      ratio = 1; // NCIS
      break;
    }
#ifdef DEBUG
    std::cout << "ChIP/Input Ratio for chr " << chrname << ": " << ratio << std::endl;
#endif
  }
  
  int32_t getbinsize() const { return binsize; }
};

class Figure {
  std::unordered_map<std::string, ChrArray> arrays;
  std::vector<SamplePairChr> pairs;
  std::vector<bed> regionBed;
  
  void loadSampleData(DROMPA::Global &p, const chrsize &chr) {
    for (auto x: p.sample) {
      arrays[x.first] = ChrArray(p, x, chr);
    }
    for (auto &x: p.samplepair) {
#ifdef DEBUG
      std::cout << "Samplepairchr " << x.argvChIP << ", " << x.argvInput << std::endl;
#endif
      pairs.emplace_back(p, x, arrays, chr.getname());
    }

#ifdef DEBUG
    std::cout << "all WigArray:" << std::endl;
    for (auto x: arrays) {
      std::cout << x.first << ", binsize " << x.second.binsize << std::endl;
    }
    std::cout << "all SamplePair:" << std::endl;
    for (auto &x: pairs) {
      std::cout << x.argvChIP << "," << x.argvInput << ", binsize " << x.getbinsize() << std::endl;
    }
#endif
  }

public:
  Figure(DROMPA::Global &p, const chrsize &chr):
    regionBed(p.drawregion.getRegionBedChr(chr.getname()))
  {}
  
  void DrawData(DROMPA::Global &p, const chrsize &chr);
  int32_t Draw(DROMPA::Global &p, const chrsize &chr) {
    if (p.drawregion.isRegionBed() && !regionBed.size()) return 0;
    loadSampleData(p, chr);
    DrawData(p, chr);
    return 1;
  }
  
};

#endif /* _DD_READFILE_H_ */
