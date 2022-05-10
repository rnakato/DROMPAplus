/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
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

    int32_t n_autosome(0);
    for (size_t i=0; i<getnchr(); ++i) {
      if (chr[i].isautosome()) {
	dist.addmp2genome(i);
	++n_autosome;
      }
    }
    if(!n_autosome) PRINTERR_AND_EXIT("No autosome appeared in the genometable file.");
  } else {
    int32_t id_longestChr = getIdLongestChr();
    genThread(dist, *this, id_longestChr, id_longestChr, prefix, sspst.isEachchr(), sspst.getNgTo());

    int32_t n_autosome(0);
    for (size_t i=0; i<getnchr(); ++i) {
      if (chr[i].isautosome()) {
	dist.addmp2genome(i);
	++n_autosome;
      }
    }
    if(!n_autosome) PRINTERR_AND_EXIT("No autosome appeared in the genometable file.");
  }

  dist.setflen(dist.name);
  dflen.setflen_ssp(dist.getnsci());

  std::cout << "\nEstimated fragment length: " << dist.getnsci() << std::endl;

  if(verbose) {
    std::string prefix2 = head + "." + typestr;
    dist.outputmpGenome(prefix2);

    setSSPstats(sspst, dist.getbackgroundUniformity(), dist.getnsc(), dist.getrlsc(), dist.getrsc());
  }

  DEBUGprint("makeProfile: " + typestr + " done.");
  DEBUGprint("strShiftProfileDROMPA done.");
  return;
}
