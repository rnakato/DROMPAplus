/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "LibraryComplexity.hpp"
#include "Mapfile.hpp"

namespace {
  using mapIntInt = std::unordered_map<int32_t, int32_t>;
  using mapStrInt = std::unordered_map<std::string, int32_t>;
  
  void hashFilterAllSingle(mapIntInt &mp, SeqStats &chr, const Strand::Strand strand, const int32_t thre)
  {
    for(auto &x: chr.getvReadref_notconst(strand)) {
      if(mp.find(x.F3) != mp.end()) {
	if(mp[x.F3] < thre) {
	  ++mp[x.F3];
	  chr.incNReadNonred(strand);
	} else {
	  x.duplicate = 1;
	  chr.incNReadRed(strand);
	}
      } else {
	mp[x.F3] = 1;
	chr.incNReadNonred(strand);
      }
    }
    return;
  }
  
  void hashFilterAllPair(mapStrInt &mp, SeqStats &chr, const Strand::Strand strand, const int32_t thre)
  {
    for(auto &x: chr.getvReadref_notconst(strand)) {
      int32_t Fmin(std::min(x.F3, x.F5));
      int32_t Fmax(std::max(x.F3, x.F5));
      std::string str(std::to_string(Fmin) + "-" + std::to_string(Fmax));
      if(mp.find(str) != mp.end()) {
	if(mp[str] < thre) {
	  ++mp[str];
	  chr.incNReadNonred(strand);
	} else {
	  x.duplicate = 1;
	  chr.incNReadRed(strand);
	}
      } else {
	mp[str] = 1;
	chr.incNReadNonred(strand);
      }
    }
    return;
  }
}

void LibComp::hashFilterCmpSingle(mapIntInt &mp, const SeqStats &chr, const Strand::Strand strand)
{
  for(auto &x: chr.getvReadref(strand)) {
    if(rand() >= r4cmp) continue;
    ++nt_all;
    if(mp.find(x.F3) != mp.end()) {
      if(mp[x.F3] < getThreshold()) {
	++mp[x.F3];
	++nt_nonred;
      } else {
	++nt_red;
      }
    } else {
      mp[x.F3] = 1;
      ++nt_nonred;
    }
  }
  return;
}

void LibComp::hashFilterCmpPair(mapStrInt &mp, const SeqStats &chr, const Strand::Strand strand)
{
  for(auto &x: chr.getvReadref(strand)) {
    if(rand() >= r4cmp) continue;
    int32_t Fmin(std::min(x.F3, x.F5));
    int32_t Fmax(std::max(x.F3, x.F5));
    std::string str(std::to_string(Fmin) + "-" + std::to_string(Fmax));
    ++nt_all;
    if(mp.find(str) != mp.end()) {
      if(mp[str] < getThreshold()) {
	++mp[str];
	++nt_nonred;
      } else {
	++nt_red;
      }
    } else {
      mp[str] = 1;
      ++nt_nonred;
    }
  }
  return;
}

void LibComp::filtering_eachchr_single(SeqStats &chr)
{
  for (auto strand: {Strand::FWD, Strand::REV}) {
    mapIntInt mp;
    hashFilterAllSingle(mp, chr, strand, getThreshold());
    mapIntInt mp2;
    hashFilterCmpSingle(mp2, chr, strand);
  }
}

void LibComp::filtering_eachchr_pair(SeqStats &chr)
{
  mapStrInt mp;
  for (auto strand: {Strand::FWD, Strand::REV}) {
    hashFilterAllPair(mp, chr, strand, getThreshold());
  }
  mapStrInt mp2;
  for (auto strand: {Strand::FWD, Strand::REV}) {
    hashFilterCmpPair(mp2, chr, strand);
  }
}

void LibComp::checkRedundantReads(SeqStatsGenome &genome)
{
  if(nofilter) {
    for (auto &x: genome.chr) x.setnread_nonread_nofilter();
    return;
  }
  
  uint64_t nread(genome.getnread(Strand::BOTH));
  setThreshold(nread, genome.getlenmpbl());
  r4cmp = getratio(ncmp, nread) * RAND_MAX;
  if(ncmp > nread) {
    std::cerr << "Warning: number of reads (" << nread << ") is < "<< static_cast<int32_t>(ncmp/NUM_1M) << " million.\n";
    lackOfRead = true;
  }

  for (auto &x: genome.chr) {
    if (ispaired) filtering_eachchr_pair(x);
    else          filtering_eachchr_single(x);
  }
  
  printf("done.\n");
  return;
}

  
