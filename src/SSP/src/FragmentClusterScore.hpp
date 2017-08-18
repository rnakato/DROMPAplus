/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _FRAGMENTCLUSTERSCORE_H_
#define _FRAGMENTCLUSTERSCORE_H_

#include <fstream>
#include <boost/bind.hpp>
#include "../common/inline.hpp"
#include "../common/BoostOptions.hpp"
#include "../common/seq.hpp"

namespace {
  const int32_t sizeOfvNeighborFrag(5000);
  //  const std::vector<int32_t> v4pnf{50, 100, 150, 500, 1000, 2000, 3000, 10000, 100000, 1000000};
  const std::vector<int32_t> v4pnf{10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 190, 200, 300, 400, 500, 700, 1000, 2000, 3000, 5000, 10000, 100000, 1000000};
}

class SeqStats;
class SeqStatsGenome;

class FCSstats {
  MyOpt::Opts opt;

  int32_t num4fcs;
  int32_t ng_from_fcs, ng_to_fcs, ng_step_fcs;
  double fcsread;
  double fcsflen;
  double fcs1k;
  double fcs10k;
  double fcs100k;

public:
  FCSstats():
    opt("Fragment cluster score",100),
    fcsread(0), fcsflen(0), fcs1k(0), fcs10k(0), fcs100k(0)
  {
    opt.add_options()
      ("ng_from_fcs",
       boost::program_options::value<int32_t>()->default_value(NUM_100K)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--ng_from_fcs")),
       "fcs start of background")
      ("ng_to_fcs",
       boost::program_options::value<int32_t>()->default_value(NUM_1M)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--ng_to_fcs")),
       "fcs end of background")
      ("ng_step_fcs",
       boost::program_options::value<int32_t>()->default_value(NUM_100K)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--ng_step_fcs")),
       "fcs step on of background")
      ;
  }
  
  void setOpts(MyOpt::Opts &allopts) {
    allopts.add(opt);
  }
  
  void setValues(const MyOpt::Variables &values) {
    DEBUGprint("FCSstats setValues...");
    
    num4fcs     = MyOpt::getVal<int32_t>(values, "num4ssp");
    ng_from_fcs = MyOpt::getVal<int32_t>(values, "ng_from_fcs");
    ng_to_fcs   = MyOpt::getVal<int32_t>(values, "ng_to_fcs");
    ng_step_fcs = MyOpt::getVal<int32_t>(values, "ng_step_fcs");
    
    DEBUGprint("FCSstats setValues done.");
  }
  void dump() {
    std::cout << boost::format("FCS background region: [%d,%d], step %d\n") % ng_from_fcs % ng_to_fcs % ng_step_fcs;
    std::cout << boost::format("Read number for FCS: %d\n") % num4fcs;
  }
  int32_t getnum4fcs()   const { return num4fcs; }
  int32_t getNgFromFCS() const { return ng_from_fcs; }
  int32_t getNgToFCS()   const { return ng_to_fcs; }
  int32_t getNgStepFCS() const { return ng_step_fcs; }
  void setfcsread(const double c) { fcsread = c; }
  void setfcsflen(const double c) { fcsflen = c; }
  void setfcs1k(const double c)   { fcs1k   = c; }
  void setfcs10k(const double c)  { fcs10k  = c; }
  void setfcs100k(const double c) { fcs100k = c; }

  void printhead(std::ofstream &out) {
    out << "FCS(read)\tFCS(flen)\tFCS(1k)\tFCS(10k)\tFCS(100k)" << std::endl;
  }
  void print(std::ofstream &out) {
    out << fcsread << "\t" << fcsflen << "\t" << fcs1k << "\t"
	<< fcs10k  << "\t" << fcs100k << std::endl;
  }
};

void makeFCSProfile(FCSstats &, const SeqStatsGenome &genome, const std::string &head, const std::string &typestr);

#endif /* _FRAGMENTCLUSTERSCORE_H_ */
