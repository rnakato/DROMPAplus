/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "WigStats.hpp"

uint32_t getWigDistThre(const std::vector<uint64_t> &wigDist, const uint64_t sum) {
  uint32_t thre(9);
  uint64_t num;
  do{
    ++thre;
    num=0;
    for(uint32_t i=0; i<thre; ++i) num += wigDist[i];
  } while(num < sum*0.8 && thre < wigDist.size() -1);
#ifdef DEBUG
  std::cout << boost::format("\nthre %1%  (%2% / %3%)\n") % thre % num % sum;
#endif
  return thre;
}

void WigStats::setWigStats(const WigArray &wigarray)
{
  double num95(wigarray.getPercentile(0.95));
  std::vector<int32_t> ar;
  int32_t wigDistSize(wigDist.size());

  for (size_t i=0; i<wigarray.size(); ++i) {
    int32_t v(wigarray[i]);
    if (v < wigDistSize) ++wigDist[v];
    if (v >= num95) continue;
    ar.emplace_back(v);
  }

  setZINBParam(ar);
}

void WigStats::peakcall(const WigArray &wigarray, const std::string chrname)
{
  if(!ave) return;
    
  int32_t ext(0);
    
  for (size_t i=0; i<wigarray.size(); ++i) {
    double val(wigarray[i]);
    double logp(getlogp(val));
      
    if (!ext) {
      if (logp > pthre) {
	vPeak.emplace_back(Peak(chrname, i, i, val, logp));
	ext=1;
      }
    } else {
      if (logp > pthre) vPeak[vPeak.size()-1].renew(i, val,logp);
      else ext=0;
    }
  }
  return;
}
