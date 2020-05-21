/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "SeqStatsDROMPA.hpp"
#include "../submodules/SSP/src/ShiftProfile.hpp"
#include "../submodules/SSP/src/ShiftProfile_p.hpp"

void SeqStatsGenome::strShiftProfile(SSPstats &sspst, const std::string &head, const bool isallchr, const bool verbose)
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
    int32_t id_longestChr = getIdLongestChr();
    genThread(dist, *this, id_longestChr, id_longestChr, prefix, sspst.isEachchr(), sspst.getNgTo());
    for (size_t i=0; i<getnchr(); ++i) {
      if (chr[i].isautosome()) dist.addmp2genome(i);
    }
  }

  dist.setflen(dist.name);
  dflen.setflen_ssp(dist.getnsci());

  if(verbose) {
    std::string prefix2 = head + "." + typestr;
    dist.outputmpGenome(prefix2);

    setSSPstats(sspst, dist.getbackgroundUniformity(), dist.getnsc(), dist.getrlsc(), dist.getrsc());
  }

  DEBUGprint("makeProfile: " + typestr + " done.");
  DEBUGprint("strShiftProfileDROMPA done.");
  return;
}
