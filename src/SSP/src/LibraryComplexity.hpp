/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _LIBRARYCOMPLEXITY_HPP_
#define _LIBRARYCOMPLEXITY_HPP_

#include <iostream>
#include <string>
#include <stdint.h>
#include <unordered_map>
#include <boost/format.hpp>
#include "SeqStats.hpp"
#include "../common/BoostOptions.hpp"

class SeqStatsGenome;

class LibComp {
  MyOpt::Opts opt;

  int32_t nofilter;
  int32_t ispaired;
  int32_t thredef;
  uint64_t ncmp;
  
  uint32_t nt_all, nt_nonred, nt_red;
  int32_t threshold;
  double r4cmp;
  bool lackOfRead;

 public:
  LibComp():
    opt("Library complexity",100),
    nt_all(0), nt_nonred(0), nt_red(0), threshold(0),
    r4cmp(0), lackOfRead(false)
  {
    opt.add_options()
      ("thre_pb",
       boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 0, "--thre_pb")),
       "PCRbias threshold (default: more than max(1 read, 10 times greater than genome average)) ")
      ("ncmp",
       boost::program_options::value<int64_t>()->default_value(NUM_10M)->notifier(boost::bind(&MyOpt::over<int64_t>, _1, 0, "--ncmp")),
       "read number for calculating library complexity")
      ("nofilter", "do not filter PCR bias")
      ;
  }

  void setValues(const MyOpt::Variables &values) {
    DEBUGprint("Libcomp setValues...");

    nofilter = values.count("nofilter");
    ispaired = values.count("pair");
    ncmp    = MyOpt::getVal<int64_t>(values, "ncmp");
    thredef = MyOpt::getVal<int32_t>(values, "thre_pb");

    DEBUGprint("Libcomp setValues done.");
  }

  void setOpts(MyOpt::Opts &allopts) {
    allopts.add(opt);
  }
  
  double getcomplexity() const { return getratio(nt_nonred, nt_all); }
  
  void setThreshold(const uint64_t nread, const uint64_t lenmpbl) {
    if(thredef) threshold = thredef;
    else {
      uint32_t thre = getratio(nread, lenmpbl) * 10;
      threshold = thre > 1 ? thre : 1;
    }
    std::cout << "Checking redundant reads: redundancy threshold " << threshold << std::endl;
  }
  int32_t getThreshold() const { return threshold; };

  void checkRedundantReads(SeqStatsGenome &genome);
  void filtering_eachchr_single(SeqStats &chr);
  void filtering_eachchr_pair(SeqStats &chr);
  void hashFilterCmpSingle(std::unordered_map<int32_t, int32_t> &mp, const SeqStats &chr, const Strand::Strand strand);
  void hashFilterCmpPair(std::unordered_map<std::string, int32_t> &mp, const SeqStats &chr, const Strand::Strand strand);

  void print(std::ofstream &out) const {
    if(lackOfRead) out << boost::format("Library complexity: (%1$.3f) (%2%/%3%)\n") % getcomplexity() % nt_nonred % nt_all;
    else           out << boost::format("Library complexity: %1$.3f (%2%/%3%)\n")   % getcomplexity() % nt_nonred % nt_all;
  }
};

#endif /* _LIBRARYCOMPLEXITY_HPP_ */
