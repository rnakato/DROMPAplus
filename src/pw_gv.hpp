/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _PW_GV_H_
#define _PW_GV_H_

#include <fstream>
#include <numeric>
#include "WigStats.hpp"
#include "GenomeCoverage.hpp"
#include "SSP/src/MThread.hpp"
#include "SSP/src/LibraryComplexity.hpp"
#include "SSP/src/Mapfile.hpp"
#include "SSP/src/ssp_shiftprofile.hpp"

namespace RPM {
  class Pnorm {
    MyOpt::Opts opt;
    std::string ntype;
    int32_t nrpm;
    double ndepth;

  public:
    Pnorm():
      opt("Total Read normalization",100)
    {
      opt.add_options()
	("ntype,n",
	 boost::program_options::value<std::string>()->default_value("NONE"),
	 "Total read normalization\n{NONE|GR|GD|CR|CD}\n   NONE: not normalize\n   GR: for whole genome, read number\n   GD: for whole genome, read depth\n   CR: for each chromosome, read number\n   CD: for each chromosome, read depth")
	("nrpm",
	 boost::program_options::value<int32_t>()->default_value(2*NUM_10M)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 0, "--nrpm")),
	 "(For GR|CR) Total read number after normalization")
	("ndepth",
	 boost::program_options::value<double>()->default_value(1.0)->notifier(boost::bind(&MyOpt::over<double>, _1, 0, "--ndepth")),
	 "(For GD|CD) Averaged read depth after normalization")
	;
    }

    void setOpts(MyOpt::Opts &allopts) {
      allopts.add(opt);
    }
    
    void setValues(const MyOpt::Variables &values) {
      DEBUGprint("Pnorm setValues...");
      ntype  = values["ntype"].as<std::string>();
      nrpm   = values["nrpm"].as<int32_t>();
      ndepth = values["ndepth"].as<double>();
      if(ntype != "NONE" && ntype != "GR" && ntype != "GD" && ntype != "CR" && ntype != "CD") PRINTERR("invalid --ntype.\n");
      DEBUGprint("Pnorm setValues done.");
    }
    void dump() const {
      std::cout << "\nTotal read normalization: " << ntype << std::endl;
      if(ntype == "GR" || ntype == "CR"){
	std::cout << "\tnormed read: " << nrpm << " for genome" << std::endl;
      }
      else if(ntype == "GD" || ntype == "CD"){
	std::cout << "\tnormed depth: " << ndepth << std::endl;
      }
      printf("\n");
    }
    
    const std::string & getType() const { return ntype; }
    int32_t getnrpm()  const { return nrpm; }
    double getndepth() const { return ndepth; }
  };
}

class Mapfile: private Uncopyable {
  int32_t on_GCnorm;
  std::string GCdir;

  bool Greekchr;

  std::string samplename;
  std::string oprefix;
  std::string obinprefix;
  std::string mpdir;
  double mpthre;

  std::vector<Peak> vPeak;
  int32_t id_longestChr;

  // GC bias
  int32_t maxGC;

 public:
  SeqStatsGenome genome;
  WigStatsGenome wsGenome;
  RPM::Pnorm rpm;
  GenomeCov::Genome gcov;

  // for SSP
  SSPstats sspst;

  class LibComp complexity;
  
  Mapfile(): Greekchr(false), mpdir(""),
    id_longestChr(0),
    maxGC(0), genome(), complexity() {}
    
  void setOpts(MyOpt::Opts &allopts) {
    genome.setOpts(allopts);
    wsGenome.setOpts(allopts);
    rpm.setOpts(allopts);
    complexity.setOpts(allopts);
    sspst.setOpts(allopts);
  }

  void setValues(const MyOpt::Variables &values);
  void dump() const {
    if(genome.isBedOn()) std::cout << "Bed file: " << genome.getbedfilename() << std::endl;
    if(mpdir != "") {
      printf("Mappability normalization:\n");
      std::cout << "\tFile directory: " << mpdir << std::endl;
      std::cout << "\tLow mappablitiy threshold: " << mpthre << std::endl;
    }
    if(on_GCnorm) {
      printf("Correcting GC bias:\n");
      std::cout << "\tChromosome directory: " << GCdir << std::endl;
    }
  }

  int32_t getIdLongestChr () const { return id_longestChr; }
  int32_t isGCnorm () const { return on_GCnorm; }
  const std::string & getSampleName() const { return samplename; }
  const std::string & getMpDir()      const { return mpdir; }

  void setmaxGC(const int32_t m) { maxGC = m; }
  int32_t getmaxGC() const {return maxGC; }

  void calcGenomeCoverage() {
    std::cout << "calculate genome coverage.." << std::flush;

    gcov.setr4cmp(genome.getnread_nonred(Strand::BOTH), genome.getnread_inbed());

    for(size_t i=0; i<genome.chr.size(); i++) {
      auto array = GenomeCov::makeGcovArray(*this, genome.chr[i], gcov.getr4cmp());
      GenomeCov::Chr chr(array, gcov.getlackOfRead());
      gcov.chr.push_back(chr);
    }
    std::cout << "done." << std::endl;
  }

  void printPeak() const {
    std::string filename = getbinprefix() + ".peak.xls";
    std::ofstream out(filename);

    vPeak[0].printHead(out);
    for(uint32_t i=0; i<vPeak.size(); ++i) {
      vPeak[i].print(out, i, wsGenome.getbinsize());
    }
  }
  const std::string & getprefix() const { return oprefix; }
  const std::string & getbinprefix() const { return obinprefix; }
  double getmpthre() const { return mpthre; }

  void addPeak(const Peak &peak) {
    vPeak.push_back(peak);
  }
  void renewPeak(const int32_t i, const double val, const double p) {
    vPeak[vPeak.size()-1].renew(i, val, p);
  }

};

#endif /* _PW_GV_H_ */
