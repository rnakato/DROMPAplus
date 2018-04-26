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
    clock_t t1,t2;
    t1 = clock();
    if(p.getSmoothing()) array.Smoothing(p.getSmoothing());
    t2 = clock();
    PrintTime(t1, t2, "Smoothing");
    t1 = clock();
    stats.setWigStats(array);
    t2 = clock();
    PrintTime(t1, t2, "WigStats");
    t1 = clock();
    stats.peakcall(array, chr.getname());
    t2 = clock();
    PrintTime(t1, t2, "peakcall");
  }
};
  
class SamplePairChr {
public:
  SamplePair &pair;
  std::vector<bed> peaks1st;
  std::vector<bed> peaks2nd;

  SamplePairChr(const DROMPA::Global &p, SamplePair &x,
		const std::unordered_map<std::string, ChrArray> &arrays,
		const std::string &chrname):
    pair(x)
  {
    setRatio(p, pair.first, arrays, chrname);
    if (x.first.peak_argv != "") peaks1st = x.first.getpeaksChr(chrname);
    if (pair.overlay) {
      setRatio(p, pair.second, arrays, chrname);
      if (x.second.peak_argv != "") peaks2nd = x.second.getpeaksChr(chrname);
    }
  }
  
  void setRatio (const DROMPA::Global &p, SamplePairParam &pair,
		 const std::unordered_map<std::string, ChrArray> &arrays, const std::string &chrname) {
    DEBUGprint("setRatio");
    if (pair.argvInput == "") return;
    pair.ratio = 1;
    
    switch (p.getNorm()) {
    case 0:
      pair.ratio = 1;
      break;
    case 1:
      pair.ratio = getratio(arrays.at(pair.argvChIP).totalreadnum, arrays.at(pair.argvInput).totalreadnum);
      break;
    case 2:
      pair.ratio = getratio(arrays.at(pair.argvChIP).totalreadnum_chr.at(chrname), arrays.at(pair.argvInput).totalreadnum_chr.at(chrname));
      break;
    case 3:
      pair.ratio = 1; // NCIS
      break;
    }
#ifdef DEBUG
    std::cout << "ChIP/Input Ratio for chr " << chrname << ": " << pair.ratio << std::endl;
#endif
  }
};

class Figure {
  std::unordered_map<std::string, ChrArray> arrays;
  std::vector<SamplePairChr> pairs;
  std::vector<bed> regionBed;
  
  void loadSampleData(DROMPA::Global &p, const chrsize &chr) {
    clock_t t1,t2;
    for (auto x: p.sample) {
      t1 = clock();
      arrays[x.first] = ChrArray(p, x, chr);
      t2 = clock();
      PrintTime(t1, t2, "ChrArray new");
    }
  
    for (auto &x: p.samplepair) {
#ifdef DEBUG
      std::cout << "Samplepairchr " << x.first.argvChIP << ", " << x.first.argvInput << std::endl;
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
      std::cout << x.pair.first.argvChIP << "," << x.pair.first.argvInput << ", binsize " << x.pair.first.getbinsize() << std::endl;
    }
#endif
  }

public:
  Figure(DROMPA::Global &p, const chrsize &chr):
    regionBed(p.drawregion.getRegionBedChr(chr.getname()))
  {}
  
  void DrawData(DROMPA::Global &p, const chrsize &chr);
  int32_t Draw(DROMPA::Global &p, const chrsize &chr) {
    clock_t t1,t2;
    if (p.drawregion.isRegionBed() && !regionBed.size()) return 0;
    t1 = clock();
    std::cout << "Loading data..";
    loadSampleData(p, chr);
    t2 = clock();
    PrintTime(t1, t2, "loadSampleData");
    t1 = clock();
    std::cout << "Drawing.." << std::endl;
    DrawData(p, chr);
    t2 = clock();
    PrintTime(t1, t2, "DrawData");
    return 1;
  }
  
};

#endif /* _DD_READFILE_H_ */
