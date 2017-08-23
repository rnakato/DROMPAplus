/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _WIGSTATS_HPP_
#define _WIGSTATS_HPP_

#include <vector>
#include <fstream>
#include <boost/bind.hpp>
#include "SSP/common/inline.hpp"
#include "SSP/common/BoostOptions.hpp"

uint32_t getWigDistThre(const std::vector<uint64_t> &, const uint64_t);

class WigArray {
  std::vector<int32_t> array;
  double geta;
  
  template <class T>
    double rmGeta(T val) const { return val/geta; } 
  template <class T>
    double addGeta(T val) const { return val*geta; }
  void checki(const size_t i) const {
    if(i>=size()) PRINTERR("Invalid i for WigArray: " << i << " > " << size());
  }

 public:
  WigArray(){}
  WigArray(const size_t num, const int32_t val): array(num, val), geta(1000.0) {}
  ~WigArray(){}

  size_t size() const { return array.size(); }
  double operator[] (const size_t i) const {
    checki(i);
    return rmGeta(array[i]);
  }

  void setval(const size_t i, const double val) {
    checki(i);
    array[i] = addGeta(val);
  }
  void addval(const size_t i, const double val) {
    checki(i);
    array[i] += addGeta(val);
  }
  void multipleval(const size_t i, const double val) {
    checki(i);
    array[i] *= val;
  }

  double getPercentile(double per) const {
    int32_t v95(MyStatistics::getPercentile(array, per));
    return rmGeta(v95);
  }

  void outputAsWig(std::ofstream &out, const int32_t binsize) const {
    for(size_t i=0; i<array.size(); ++i) {
      if(array[i]) out << boost::format("%1%\t%2%\n") % (i*binsize+1) % rmGeta(array[i]);
    }
  }
  void outputAsBedGraph(std::ofstream &out, const int32_t binsize, const std::string &name, const uint64_t chrend) const {
    uint64_t e;
    for(size_t i=0; i<array.size(); ++i) {
      if(i==array.size() -1) e = chrend; else e = (i+1) * binsize;
      if(array[i]) out << boost::format("%1% %2% %3% %4%\n") % name % (i*binsize) % e % rmGeta(array[i]);
    }
  }
  void outputAsBinary(std::ofstream &out) const {
    for(size_t i=0; i<array.size(); ++i) out.write((char *)&array[i], sizeof(int32_t));
  }
  void readBinary(std::ifstream &in, const int32_t nbin) const {
    for(int32_t i=0; i<nbin; ++i) in.read((char *)&array[i], sizeof(int32_t));
  }
  void dump() const {
    for (auto x: array) std::cout << x << std::endl;
  }
};

class WigStats {
  enum{ n_wigDist=200 };

  uint64_t len;
  int32_t binsize;
  uint64_t sum;

 public:
  double ave, var, nb_p, nb_n, nb_p0;
  std::vector<uint64_t> wigDist;

  WigStats(uint64_t l, int32_t b): len(l), binsize(b), sum(0), ave(0), var(0), nb_p(0), nb_n(0), nb_p0(0),
    wigDist(n_wigDist,0) {}

  uint64_t getsum() const { return sum; };
  int32_t getnbin() const { return binsize ? len/binsize +1: 0; }
  double getPoisson(const int32_t i) const {
    if(ave) return _getPoisson(i, ave);
    else return 0;
  }
  double getNegativeBinomial(const int32_t i) const {
    return _getNegativeBinomial(i, nb_p, nb_n);
  }
  double getZINB(const int32_t i) const {
    if(ave) return _getZINB(i, nb_p, nb_n, nb_p0);
    else return 0;
  }
  void getWigStats(const WigArray &wigarray) {
    double num95 = wigarray.getPercentile(0.95);
    
    int32_t size = wigDist.size();
    std::vector<int32_t> ar;
    for(size_t i=0; i<wigarray.size(); ++i) {
      //int32_t v(wigarray.getval(i));
      int32_t v(wigarray[i]);
      if(v<0) std::cout << sum << "xxx" << v << std::endl;
      ++sum;
      if(v < size) ++wigDist[v];
      if(v >= num95) continue;
      ar.push_back(v);
    }

    MyStatistics::moment<int32_t> mm(ar, 0);
    ave = mm.getmean();
    var = mm.getvar();
    nb_p = var ? ave/var : 0;
    if(nb_p>=1) nb_p = 0.9;
    if(nb_p<=0) nb_p = 0.1; 
    nb_n = ave * nb_p /(1 - nb_p);

    //    std::cout << ave << "\t" << var << "\t" << nb_p << "\t" << nb_n<< std::endl;
    if(ave) estimateParam();
  }

  double getpWig(const int32_t i) const { return getratio(wigDist[i], sum); }

  void estimateParam() {
    uint32_t thre = getWigDistThre(wigDist, sum);
    double par[thre+1];
    par[0] = thre;
    for(uint32_t i=0; i<thre; ++i) par[i+1] = getpWig(i);
    iterateZINB(&par, nb_p, nb_n, nb_p, nb_n, nb_p0);
  }
  void printwigDist(std::ofstream &out, const int32_t i) const {
    out << boost::format("%1%\t%2%\t") % wigDist[i] % getpWig(i);
  }
  void printPoispar(std::ofstream &out) const {
    out << boost::format("%1$.3f\t%2$.3f\t") % ave % var;
  }
  void printZINBpar(std::ofstream &out) const {
    out << boost::format("%1%\t%2%\t%3%") % nb_p % nb_n % nb_p0;
  }
};

enum class WigType {
  BINARY,
  COMPRESSWIG,
  UNCOMPRESSWIG,
  BEDGRAPH,
  BIGWIG,
  NONE,
  WIGTYPENUM
};

class WigStatsGenome {
  enum{ n_wigDist=200 };
  std::vector<std::string> strType = {"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};

  MyOpt::Opts opt;
  int32_t binsize;
  int32_t rcenter;
  WigType type;

public:
  std::vector<WigStats> chr;
  std::vector<uint64_t> wigDist;
  double ave, var, nb_p, nb_n, nb_p0;

  WigStatsGenome():
    opt("Wigarray",100),
    wigDist(n_wigDist,0),
    ave(0), var(0), nb_p(0), nb_n(0), nb_p0(0)
  {
  opt.add_options()
    ("wigformat",
     boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::range<int32_t>, _1, 0, static_cast<int>(WigType::WIGTYPENUM) -2, "--wigformat")),
     "Output format\n   0: binary (.bin)\n   1: compressed wig (.wig.gz)\n   2: uncompressed wig (.wig)\n   3: bedGraph (.bedGraph)\n   4: bigWig (.bw)")
    ("binsize,b",
     boost::program_options::value<int32_t>()->default_value(50)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--binsize")),
     "bin size")
    ("rcenter",
     boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 0, "--rcemter")),
     "consider length around the center of fragment")
    ;
  }

  void setOpts(MyOpt::Opts &allopts) {
    allopts.add(opt);
  }
  void setValues(const MyOpt::Variables &values) {
    binsize = MyOpt::getVal<int32_t>(values, "binsize");
    rcenter = MyOpt::getVal<int32_t>(values, "rcenter");
    type    = static_cast<WigType>(MyOpt::getVal<int32_t>(values, "wigformat"));
  }
  void dump() const {
    std::cout << "\tOutput format: " << strType[static_cast<int32_t>(type)] << std::endl;
    std::cout << "Binsize: " << binsize << " bp" << std::endl;
  }

  int32_t getbinsize() const { return binsize; }
  int32_t getrcenter() const { return rcenter; }
  WigType getWigType() const { return type; }
  
  uint64_t getsum() const {
    uint64_t sum(0);
    for(auto &x: chr) sum += x.getsum();
    return sum;
  }
  int32_t getnbin() const {
    int32_t nbin(0);
    for(auto &x: chr) nbin += x.getnbin();
    return nbin;
  }
  
  void addWigArray(const int32_t id, const WigArray &array) {
    chr[id].getWigStats(array);
    for(size_t i=0; i<wigDist.size(); ++i) {
      wigDist[i] += chr[id].wigDist[i];
    }
  }
  double getpWig(const int32_t i) const { return getratio(wigDist[i], getsum());}

  double getZINB(const int32_t i) const {
    if(ave) return _getZINB(i, nb_p, nb_n, nb_p0);
    else return 0;
  }
  void printwigDist(std::ofstream &out, const int32_t i) const {
    out << boost::format("%1%\t%2%\t") % wigDist[i] % getpWig(i);
  }
  void printPoispar(std::ofstream &out) const {
    out << boost::format("%1$.3f\t%2$.3f\t") % ave % var;
  }
  void printZINBpar(std::ofstream &out) const {
    out << boost::format("%1%\t%2%\t%3%") % nb_p % nb_n % nb_p0;
  }

  void estimateZINB(const int32_t id) {
    uint32_t thre = getWigDistThre(wigDist, getsum());
    double par[thre+1];
    par[0] = thre;
    for(uint32_t i=0; i<thre; ++i) par[i+1] = getpWig(i);

    iterateZINB(&par, chr[id].nb_p, chr[id].nb_n, nb_p, nb_n, nb_p0);

    return;
  }
};

#endif /* _WIGSTATS_HPP_ */
