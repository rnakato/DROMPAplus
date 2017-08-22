/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_DRAW_H_
#define _DD_DRAW_H_

#include "dd_class.hpp"
#include "dd_readfile.hpp"

class ChrArray {
public:
  int32_t binsize;
  int32_t nbin;
  WigArray array;
    
  ChrArray(){}
  /*    Array(const std::string &filename, const int32_t b, const WigType &iftype, const chrsize &chr):
	binsize(b), nbin(chr.getlen()/binsize +1),
	array(readInputData(filename, binsize, nbin, iftype, chr))
	{}*/
  void setValue(const std::string &filename, const int32_t b, const WigType &iftype, const chrsize &chr) {
    binsize = b;
    nbin = chr.getlen()/binsize +1;
    array = readInputData(filename, binsize, nbin, iftype, chr);
  }
};
  
class SamplePairChr {
  int32_t binsize;
public:
  std::string argvChIP;
  std::string argvInput;
  std::string label;

  SamplePairChr(const std::string &argvChIP, const std::string &argvInput, const std::string &label, const int32_t &b):
    binsize(b), argvChIP(argvChIP), argvInput(argvInput), label(label) {}
  int32_t getbinsize() const { return binsize; }
};

class Figure {
  std::unordered_map<std::string, ChrArray> arrays;
  std::vector<SamplePairChr> pairs;
  std::vector<bed> regionBed;
  
  void loadSampleData(DROMPA::Global &p, const chrsize &chr) {
    for(auto x: p.sample) {
      arrays[x.first] = ChrArray();
      arrays[x.first].setValue(x.first, x.second.getbinsize(), x.second.getiftype(), chr);
    }
    for(auto itr = p.samplepair.begin(); itr != p.samplepair.end(); ++itr) {
      pairs.emplace_back(itr->argvChIP, itr->argvInput, itr->label, arrays[itr->argvChIP].binsize);
    }

#ifdef DEBUG
    std::cout << "all WigArray:" << std::endl;
    for(auto x: arrays) {
      std::cout << x.first << ", binsize " << x.second.binsize << std::endl;
    }
    std::cout << "all SamplePair:" << std::endl;
    for(auto itr = pairs.begin(); itr != pairs.end(); ++itr) {
      std::cout << itr->argvChIP << "," << itr->argvInput << ", binsize " << itr->getbinsize() << std::endl;
    }
#endif
  }

public:
  Figure(DROMPA::Global &p, const chrsize &chr):
    regionBed(p.drawregion.getRegionBedChr(chr.getname()))
  {}
  
  void DrawData(DROMPA::Global &p, const chrsize &chr);
  int32_t Draw(DROMPA::Global &p, const chrsize &chr) {
    if(p.drawregion.isRegionBed() && !regionBed.size()) return 0;
    loadSampleData(p, chr);
    DrawData(p, chr);
    return 1;
  }
  
};


#endif /* _DD_READFILE_H_ */
