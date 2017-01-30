/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _PW_GV_H_
#define _PW_GV_H_

#include <fstream>
#include <numeric>
#include <boost/thread.hpp>
#include "WigStats.hpp"
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
  int32_t on_bed;
  int32_t on_GCnorm;
  std::string bedfilename;
  std::string GCdir;

  bool Greekchr;

  std::string samplename;
  std::string oprefix;
  std::string obinprefix;
  std::string mpdir;
  double mpthre;

  bool lackOfRead4GenomeCov;
  bool lackOfRead4FragmentVar;
  std::vector<Peak> vPeak;
  int32_t id_longestChr;

  // GC bias
  int32_t maxGC;

 public:
  SeqStatsGenome genome;
  WigStatsGenome wsGenome;
  RPM::Pnorm rpm;

  // for SSP
  SSPstats sspst;

  class LibComp complexity;
  
  Mapfile(): Greekchr(false), mpdir(""),
    lackOfRead4GenomeCov(false),
    lackOfRead4FragmentVar(false),
    id_longestChr(0),
    maxGC(0), genome(), complexity() {}
    
  void setOpts(MyOpt::Opts &allopts) {
    genome.setOpts(allopts);
    wsGenome.setOpts(allopts);
    rpm.setOpts(allopts);
    complexity.setOpts(allopts);
    sspst.setOpts(allopts);
  }
  void setValues(const MyOpt::Variables &values) {
    DEBUGprint("Mapfile setValues...");
    
    on_bed = values.count("bed");
    if(on_bed) bedfilename = values["bed"].as<std::string>();
    on_GCnorm = values.count("genome");
    if(on_GCnorm) GCdir = values["genome"].as<std::string>();
    genome.setValues(values);
    wsGenome.setValues(values);

    for(auto itr = genome.chr.begin(); itr != genome.chr.end(); ++itr) {
      wsGenome.chr.push_back(WigStats(itr->getlen(), wsGenome.getbinsize()));
    }
    
    rpm.setValues(values);
    complexity.setValues(values);
    sspst.setValues(values);

    samplename = values["output"].as<std::string>();
    id_longestChr = setIdLongestChr(genome);
    oprefix = values["odir"].as<std::string>() + "/" + values["output"].as<std::string>();
    obinprefix = oprefix + "." + IntToString(values["binsize"].as<int32_t>());

    if (values.count("mp")) mpdir = values["mp"].as<std::string>();
    mpthre = values["mpthre"].as<double>();
    
    DEBUGprint("Mapfile setValues done.");
  }

  void dump() const {
    if(on_bed) std::cout << "Bed file: " << bedfilename << std::endl;
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
  int32_t isBedOn  () const { return on_bed; }
  const std::string & getSampleName() const { return samplename; }
  const std::string & getMpDir()   const { return mpdir; }

  void setmaxGC(const int32_t m) { maxGC = m; }
  int32_t getmaxGC() const {return maxGC; }

  void lackOfRead4GenomeCov_on() { lackOfRead4GenomeCov = true; }
  bool islackOfRead4GenomeCov() const { return lackOfRead4GenomeCov; };
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
