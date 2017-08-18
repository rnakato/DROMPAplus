/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _MTHREAD_HPP_
#define _MTHREAD_HPP_

#include <cstdint>
#include <vector>

namespace MyMthread {
  class chrrange {
  public:
    uint32_t s;
    uint32_t e;
  chrrange(uint32_t start, uint32_t end): s(start), e(end) {}
  };

  template <class T>
  std::vector<chrrange> getVsepchr(const uint64_t genomelen, const std::vector<T> &chr, const int32_t numthreads)
  {
    std::vector<chrrange> vsep;
    
    uint32_t sepsize = genomelen/numthreads;
    for(uint32_t i=0; i<chr.size(); ++i) {
      uint32_t s = i;
      uint64_t l(0);
      while(l < sepsize && i<chr.size()) {
	l += chr[i].getlen();
	i++;
      }
      i--;
      uint32_t e = i;
      vsep.emplace_back(s, e);
    }
    return vsep;
  }
}

#endif /* _MTHREAD_HPP_ */
