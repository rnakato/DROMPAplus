/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_READFILE_H_
#define _DD_READFILE_H_

#include <unordered_map>
#include <string>
#include <vector>
#include "dd_gv.hpp"
#include "../submodules/SSP/common/seq.hpp"

WigArray loadWigData(const std::string &filename, const SampleInfo &x, const chrsize &, const int32_t binnum_for_localpoisson);

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
    binsize(x.second.getbinsize()), 
    nbin(chr.getlen()/binsize +1),
    array(loadWigData(x.first, x.second, chr, p.thre.binnum_for_localpoisson)),
    stats(nbin, binsize),
    totalreadnum(x.second.gettotalreadnum()),
    totalreadnum_chr(x.second.gettotalreadnum_chr())
  {
    clock_t t1,t2;
    t1 = clock();
    if(p.getSmoothing()) array.Smoothing(p.getSmoothing(), p.getSmoothingtype());
    t2 = clock();
    PrintTime(t1, t2, "Smoothing");
    t1 = clock();
    stats.setWigStats(array);
    t2 = clock();
    PrintTime(t1, t2, "WigStats");
  }
};

class vChrArray {
  const chrsize &chr;
  std::unordered_map<std::string, ChrArray> arrays;

public:
  vChrArray(const DROMPA::Global &p, const chrsize &_chr);

  const ChrArray & getArray(const std::string &str) const {
    return arrays.at(str);
  }
  const chrsize & getchr() const { return chr; }
  int32_t getchrlen() const { return chr.getlen(); }
};

#endif /* _DD_READFILE_H_ */
