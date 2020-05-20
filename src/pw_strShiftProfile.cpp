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
  for(size_t i=0; i<genome.chr.size(); ++i) {
    if (lenmax < genome.chr[i].getlenmpbl()) {
      lenmax = genome.chr[i].getlenmpbl();
      id = i;
    }
  }
  return id;
}

void makeProfile_forDROMPA(SSPstats &sspst, SeqStatsGenome &genome, const std::string &head, const std::string &typestr)
{
  DEBUGprint("makeProfile: " + typestr);
  shiftJacBit dist(sspst, genome);
  dist.printStartMessage();

  int32_t id_longestChr = setIdLongestChr(genome);

  std::string prefix(head + "." + typestr);
  genThread(dist, genome, id_longestChr, id_longestChr, prefix, sspst.isEachchr(), sspst.getNgTo());

  for (size_t i=0; i<genome.chr.size(); ++i) {
    if (genome.chr[i].isautosome()) dist.addmp2genome(i);
  }

  dist.setflen(dist.name);
  genome.dflen.setflen_ssp(dist.getnsci());

  if(sspst.getNgTo() < 0) return;

  std::string prefix2 = head + "." + typestr;
  dist.outputmpGenome(prefix2);

  if (typestr == "jaccard") setSSPstats(sspst, dist.getbackgroundUniformity(), dist.getnsc(), dist.getrlsc(), dist.getrsc());

  DEBUGprint("makeProfile: " + typestr + " done.");
  return;
}
void strShiftProfile(SSPstats &sspst, SeqStatsGenome &genome, const std::string &head)
{
  DEBUGprint("strShiftProfileDROMPA...");
//  if (genome.dflen.isallchr()) makeProfile<shiftJacBit>(sspst, genome, head, "jaccard");
//  else makeProfile_forDROMPA(sspst, genome, head, "jaccard");

 // use the longest chromosome only
  makeProfile_forDROMPA(sspst, genome, head, "jaccard");

  DEBUGprint("strShiftProfileDROMPA done.");
  return;
}
