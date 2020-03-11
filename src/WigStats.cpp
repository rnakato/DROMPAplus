/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "WigStats.hpp"

double getlogp_Poisson(const double val, const double myu);

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
  //  double min(wigarray.getMinValue());
  std::vector<int32_t> ar;
  int32_t wigDistSize(wigDist.size());

  //  std::cout << "num95  "<< num95 << "  min   " << min << std::endl;

  // 負の値は暫定的に無視
  for (size_t i=0; i<wigarray.size(); ++i) {
    int32_t v(wigarray[i]);
    if (v >= num95 || v < 0) continue;
    if (v < wigDistSize) ++wigDist[v];
    ar.emplace_back(v);
  }

  //  setZINBParam(ar);
}

/*void WigStats::peakcall_onlyChIP(const WigArray &wigarray, const std::string &chrname)
{
  int32_t ext(0);
  for (size_t i=0; i<wigarray.size(); ++i) {
    double val(wigarray[i]);
    double logp_inter(getlogp_Poisson(val, wigarray.getLocalAverage(i, binsize)));

    if (!ext) {
      if (logp_inter > pthre) {
	vPeak.emplace_back(Peak(chrname, i, i, val, logp_inter));
	ext=1;
      }
    } else {
      if (logp_inter > pthre) vPeak[vPeak.size()-1].renew(i, val, logp_inter);
      else ext=0;
    }
  }
  return;
}
*/
