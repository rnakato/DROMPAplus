/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _GCNORMALIZATION_HPP_
#define _GCNORMALIZATION_HPP_

#include <numeric>
#include <boost/bind.hpp>
#include "../submodules/SSP/common/BoostOptions.hpp"
#include "../submodules/SSP/common/util.hpp"
#include "../submodules/SSP/common/inline.hpp"

class bed;
class SeqStats;
class SeqStatsGenome;

class GCnorm {
  MyOpt::Opts opt;
  
  int32_t on_GCnorm;
  std::string GCdir;
  int32_t flen4gc;
  int32_t gcdepthoff;

public:
  GCnorm():
    opt("GC normalization",100),
    on_GCnorm(0), GCdir(""),
    flen4gc(0), gcdepthoff(0)
  {
    opt.add_options()
      ("chrdir", boost::program_options::value<std::string>(),
       "Chromosome directory of reference genome sequence for GC content estimation")
      ("flen4gc",
       boost::program_options::value<int32_t>()->default_value(120)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 0, "--flen4gc")),
       "Fragment length for calculation of GC distribution")
      ("gcdepthoff", "ignore to consider depth of GC contents")
      ;
  }
  
  void setOpts(MyOpt::Opts &allopts) {
    allopts.add(opt);
  }
  void setValues(const MyOpt::Variables &values) {
    on_GCnorm  = values.count("chrdir");
    if(on_GCnorm) GCdir = MyOpt::getVal<std::string>(values, "chrdir");
    flen4gc    = MyOpt::getVal<int32_t>(values, "flen4gc");
    gcdepthoff = values.count("gcdepthoff");
  }
  
  const std::string & getGCdir() const { return GCdir; }
  int32_t isGcNormOn()   const { return on_GCnorm; }
  int32_t getflen4gc()   const { return flen4gc; }
  int32_t isGcDepthOff() const { return gcdepthoff; }

};


class GCdist {
  int32_t flen;
  int32_t flen4gc;

  std::vector<int32_t> DistGenome;
  std::vector<int32_t> DistRead;
  std::vector<double> GCweight;

  double getPropGenome(const int32_t i) {
    double GsumGC = accumulate(DistGenome.begin(), DistGenome.end(), 0);
    return getratio(DistGenome[i], GsumGC);
  }
  double getPropRead(const int32_t i) {
    double RsumGC = accumulate(DistRead.begin(), DistRead.end(), 0);
    return getratio(DistRead[i], RsumGC);
  }
  double getPropDepth(const int32_t i) {
    return getratio(DistRead[i], DistGenome[i]);
  }
  void makeGCweightDist(const int32_t);

public:
  GCdist(const int32_t l, GCnorm &gc);
  void calcGCdist(const SeqStats &chr, const GCnorm &gc, const std::string &mpdir, const int32_t isBedOn, const std::vector<bed> &vbed);

  int32_t getmaxGC() const { return getmaxi(DistRead); }
  double getGCweight(const int32_t i) const { return GCweight[i]; }
  void outputGCweightDist(const std::string &filename);
  
  int32_t getflen() const { return flen; }
  int32_t getflen4gc() const { return flen4gc; }
};

void weightRead(SeqStatsGenome &, GCdist &, const std::string &);


#endif /* _GCNORMALIZATION_HPP_ */
