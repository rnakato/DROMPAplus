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
