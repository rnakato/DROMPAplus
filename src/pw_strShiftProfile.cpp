/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "pw_strShiftProfile.hpp"
#include "SeqStatsDROMPA.hpp"
#include "../submodules/SSP/src/ShiftProfile.hpp"
#include "../submodules/SSP/src/ShiftProfile_p.hpp"

int32_t setIdLongestChr(const SeqStatsGenome &genome)
{
  int32_t id(0);
  uint64_t lenmax(0);
  for(size_t i=0; i<genome.getnchr(); ++i) {
    if (lenmax < genome.chr[i].getlenmpbl()) {
      lenmax = genome.chr[i].getlenmpbl();
      id = i;
    }
  }
  return id;
}

void SeqStatsGenome::strShiftProfile(SSPstats &sspst, const std::string &head, const bool isallchr)
{
  DEBUGprint("strShiftProfileDROMPA...");

  std::string typestr("jaccard");

  DEBUGprint("makeProfile: " + typestr);
  shiftJacBit dist(sspst, *this);
  dist.printStartMessage();

  std::string prefix(head + "." + typestr);

  if (isallchr) {
    boost::thread_group agroup;
    for (size_t i=0; i<vsepchr.size(); ++i) {
      agroup.create_thread(bind(genThread<shiftJacBit>, boost::ref(dist), boost::cref(*this), vsepchr[i].s, vsepchr[i].e, boost::cref(prefix), sspst.isEachchr(), sspst.getNgTo()));
    }
    agroup.join_all();

    for (size_t i=0; i<getnchr(); ++i) {
      if (chr[i].isautosome()) dist.addmp2genome(i);
    }
  } else {
    int32_t id_longestChr = setIdLongestChr(*this);
    genThread(dist, *this, id_longestChr, id_longestChr, prefix, sspst.isEachchr(), sspst.getNgTo());
    for (size_t i=0; i<getnchr(); ++i) {
      if (chr[i].isautosome()) dist.addmp2genome(i);
    }
  }

  dist.setflen(dist.name);
  dflen.setflen_ssp(dist.getnsci());

  if(sspst.getNgTo() < 0) return;

  std::string prefix2 = head + "." + typestr;
  dist.outputmpGenome(prefix2);

  setSSPstats(sspst, dist.getbackgroundUniformity(), dist.getnsc(), dist.getrlsc(), dist.getrsc());

  DEBUGprint("makeProfile: " + typestr + " done.");


  DEBUGprint("strShiftProfileDROMPA done.");
  return;
}
