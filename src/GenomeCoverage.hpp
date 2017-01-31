/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _GENOMECOVERAGE_HPP_
#define _GENOMECOVERAGE_HPP_

#include <iostream>
#include <fstream>
#include <stdint.h>
#include <boost/format.hpp>
#include "SSP/common/inline.hpp"
#include "SSP/src/BpStatus.hpp"

class SeqStats;
class Mapfile;

namespace GenomeCov {
  std::vector<BpStatus> makeGcovArray(const Mapfile &, const SeqStats &chr, const double r4cmp);

  class gvStats {
    virtual uint64_t getnbp() const = 0;
    virtual uint64_t getncov() const = 0;
    virtual uint64_t getncovnorm() const = 0;
    
  protected:
    bool lackOfRead;
    
  public:
    gvStats(const bool b): lackOfRead(b) {}
    virtual ~gvStats(){}

    void printstats(std::ofstream &out) const {
      if(lackOfRead) {
	std::cout<< getncov() << getncovnorm() << getnbp() << std::endl;
 	out << boost::format("%1$.3f\t(%2$.3f)\t")
	  % getratio(getncov(),     getnbp())
	  % getratio(getncovnorm(), getnbp());
      } else {
	out << boost::format("%1$.3f\t%2$.3f\t")
	  % getratio(getncov(),     getnbp())
	  % getratio(getncovnorm(), getnbp());
      }
    }
  };
  
  class Chr: public gvStats {
    uint64_t nbp, ncov, ncovnorm;
    
  public:
    Chr(const std::vector<BpStatus> &array, const bool b):
      gvStats(b), nbp(0), ncov(0), ncovnorm(0)
    {
      for(auto &x: array) {
	if(x >= BpStatus::MAPPABLE)     ++nbp;      // MAPPABLE || COVREAD_ALL || COVREAD_NORM
	if(x >= BpStatus::COVREAD_ALL)  ++ncov;     // COVREAD_ALL || COVREAD_NORM
	if(x == BpStatus::COVREAD_NORM) ++ncovnorm;
      }
    }
    
    uint64_t getnbp()       const { return nbp; }
    uint64_t getncov()      const { return ncov; }
    uint64_t getncovnorm()  const { return ncovnorm; }
  };

  class Genome: public gvStats {
    enum { numGcov=5000000 };
    double r4cmp;
    
  public:
    std::vector<Chr> chr;
    Genome(): gvStats(false), r4cmp(0) {}
    
    void setr4cmp (const uint64_t nread_nonred, const uint64_t nread_inbed) {
      // ignore peak region
      double r = getratio(numGcov, nread_nonred - nread_inbed);
      if(r>1){
	std::cerr << "Warning: number of reads is < "<< numGcov << " for GenomeCoverage. \n";
	lackOfRead = true;
      }
      r4cmp = r*RAND_MAX;
    }
    double getr4cmp() const { return r4cmp; }
    bool getlackOfRead() const { return lackOfRead; }

    uint64_t getnbp() const {
      uint64_t nbp(0);
      for(auto &x: chr) nbp += x.getnbp();
      return nbp;
    }
    uint64_t getncov() const {
      uint64_t ncov(0);
      for(auto &x: chr) ncov += x.getncov();
      return ncov;
    }
    uint64_t getncovnorm() const {
      uint64_t ncovnorm(0);
      for(auto &x: chr) ncovnorm += x.getncovnorm();
      return ncovnorm;
    }

    void printhead(std::ofstream &out) const {
      out << "gcov (Raw)\tgcov (for " << static_cast<int>(numGcov/NUM_1M) << "M reads)\t";
    }
  };
}

#endif /* _GENOMECOVERAGE_HPP_ */
