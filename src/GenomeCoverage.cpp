/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <algorithm>
#include "GenomeCoverage.hpp"
#include "ReadBpStatus.hpp"
#include "pw_gv.hpp"
#include "SSP/src/SeqStats.hpp"

namespace GenomeCov {
  std::vector<BpStatus> makeGcovArray(const Mapfile &p, const SeqStats &chr, const double r4cmp)
  {
    std::vector<BpStatus> array;
    if(p.getMpDir() != "") array = readMpbl_binary(p.getMpDir(), ("chr" + chr.getname()), chr.getlen());
    else array = readMpbl_binary(chr.getlen());
    if(p.isBedOn()) arraySetBed(array, chr.getname(), p.genome.getvbedref());

    for (auto strand: {Strand::FWD, Strand::REV}) {
      for (auto &x: chr.getvReadref(strand)) {
	if(x.duplicate) continue;
      
	BpStatus val;
	if(rand() >= r4cmp) val = BpStatus::COVREAD_ALL;
	else                val = BpStatus::COVREAD_NORM;

	int32_t s(std::max(0, std::min(x.F3, x.F5)));
	int32_t e(std::min(std::max(x.F3, x.F5), static_cast<int32_t>(chr.getlen()-1)));
	std::cout << static_cast<int>(val) << "\t"<< x.F3<< "\t"<< x.F5<< "\t"<< s<< "\t"<< e<<std::endl;
	if(s >= chr.getlen() || e < 0) {
	  std::cerr << "Warning: " << chr.getname() << " read " << s <<"-"<< e << " > array size " << chr.getlen() << std::endl;
	}
	for(int32_t i=s; i<=e; ++i) {
	  if(array[i]==BpStatus::MAPPABLE) {
	    if(rand() >= r4cmp) array[i] = val;
	  }
	}
      }
    }
    return array;
  }

}
