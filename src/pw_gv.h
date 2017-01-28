/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _PW_GV_H_
#define _PW_GV_H_

#include <fstream>
#include <numeric>
#include <boost/thread.hpp>
#include "SSP/src/readdata.h"
#include "SSP/src/mthread.h"
#include "SSP/src/LibraryComplexity.hpp"
#include "SSP/src/Mapfile.hpp"
#include "SSP/src/ssp_shiftprofile.h"

//namespace PW {
class Mapfile: private Uncopyable {
  bool Greekchr;

  std::string samplename;
  std::string oprefix;
  std::string obinprefix;

  bool lackOfRead4GenomeCov;
  bool lackOfRead4FragmentVar;
  std::vector<Peak> vPeak;

  // GC bias
  int32_t maxGC;

 public:
  SeqStatsGenome genome;
  WigStats wsGenome;
  std::vector<SeqWigStats>::iterator lchr; // longest chromosome

  // for SSP
  SSPstats sspst;

  class LibComp complexity;
  // Wigdist
  int32_t nwigdist;
  std::vector<int32_t> wigDist;
  
 Mapfile(): Greekchr(false),
    lackOfRead4GenomeCov(false),
    lackOfRead4FragmentVar(false),
    maxGC(0), genome(), complexity() {}
    
  void setOpts(MyOpt::Opts &allopts) {
    genome.setOpts(allopts);
    complexity.setOpts(allopts);
    sspst.setOpts(allopts);
  }
  void setValues(const MyOpt::Variables &values) {
    genome.setValues(values);
    complexity.setValues(values);
    samplename = values["output"].as<std::string>();
    lchr = setlchr(genome);
    oprefix = values["odir"].as<std::string>() + "/" + values["output"].as<std::string>();
    obinprefix = oprefix + "." + IntToString(values["binsize"].as<int32_t>());
  }

  void setmaxGC(const int32_t m) { maxGC = m; }
  int32_t getmaxGC() const {return maxGC; }

  void lackOfRead4GenomeCov_on() { lackOfRead4GenomeCov = true; }
  bool islackOfRead4GenomeCov() const { return lackOfRead4GenomeCov; };
  void printPeak(const MyOpt::Variables &values) const {
    std::string filename = getbinprefix() + ".peak.xls";
    std::ofstream out(filename);

    vPeak[0].printHead(out);
    for(uint32_t i=0; i<vPeak.size(); ++i) {
      vPeak[i].print(out, i, values["binsize"].as<int32_t>());
    }
  }
  std::string getprefix() const { return oprefix; }
  std::string getbinprefix() const { return obinprefix; }

  void addPeak(const Peak &peak) {
    vPeak.push_back(peak);
  }
  void renewPeak(const int32_t i, const double val, const double p) {
    vPeak[vPeak.size()-1].renew(i, val, p);
  }

  void estimateZINB() {
    int32_t thre = wsGenome.getwigDistthre();
    double par[thre+1];
    par[0] = thre;
    for(int32_t i=0; i<thre; ++i) par[i+1] = wsGenome.pwigDist[i];

    iterateZINB(&par, lchr->ws.nb_p, lchr->ws.nb_n, wsGenome.nb_p, wsGenome.nb_n, wsGenome.nb_p0);

    return;
  }
};
//}

#endif /* _PW_GV_H_ */
