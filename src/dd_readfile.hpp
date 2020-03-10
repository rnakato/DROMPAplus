/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_READFILE_H_
#define _DD_READFILE_H_

#include <unordered_map>
#include <string>
#include <vector>
#include "dd_gv.hpp"
#include "../submodules/SSP/common/seq.hpp"

WigArray loadWigData(const std::string &filename, const SampleInfo &x, const chrsize &chr);

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
	   const std::pair<const std::string, SampleInfo> &x,
	   const chrsize &chr):
    binsize(x.second.getbinsize()), nbin(chr.getlen()/binsize +1),
    array(loadWigData(x.first, x.second, chr)),
    stats(nbin, binsize, p.thre.pthre_inter),
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

using ChrArrayMap = std::unordered_map<std::string, ChrArray>;

class vChrArray {
  const chrsize &chr;
  ChrArrayMap arrays;

  void loadvChrArray(const DROMPA::Global &p) {
    std::cout << "Load sample data..";
    for (auto &x: p.vsinfo.getarray()) {
      clock_t t1,t2;
      t1 = clock();
      arrays[x.first] = ChrArray(p, x, chr);
      t2 = clock();
      PrintTime(t1, t2, "ChrArray new");
    }
  }

public:
  vChrArray(const DROMPA::Global &p, const chrsize &_chr):
    chr(_chr)
  {
    loadvChrArray(p);

#ifdef DEBUG
    std::cout << "all WigArray:" << std::endl;
    for (auto &x: arrays) {
      std::cout << x.first << ", binsize " << x.second.binsize << std::endl;
    }
    std::cout << "all SamplePair:" << std::endl;
    for (auto &x: p.samplepair) {
	std::cout << x.first.argvChIP << "," << x.first.argvInput
		  << ", binsize " << x.first.getbinsize() << std::endl;
    }
#endif
  }

  const ChrArray & getArray(const std::string &str) const {
    return arrays.at(str);
  }
  const chrsize & getchr() const { return chr; }
  int32_t getchrlen() const { return chr.getlen(); }
};


#endif /* _DD_READFILE_H_ */
