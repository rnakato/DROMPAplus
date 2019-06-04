/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _WIGSTATS_HPP_
#define _WIGSTATS_HPP_

#include <vector>
#include <fstream>
#include <boost/bind.hpp>
#include "../submodules/SSP/common/util.hpp"
#include "../submodules/SSP/common/inline.hpp"
#include "../submodules/SSP/common/BoostOptions.hpp"
#include "../submodules/SSP/common/BedFormat.hpp"
#include "../submodules/SSP/src/SeqStats.hpp"

uint32_t getWigDistThre(const std::vector<uint64_t> &, const uint64_t);

enum class WigType {
  //  BINARY,
  COMPRESSWIG,
  UNCOMPRESSWIG,
  BEDGRAPH,
  BIGWIG,
  NONE,
  WIGTYPENUM
};

class WigArray;

class WigStats {
  enum {WIGDISTNUM=200};

  double ave, var, nb_p0;
  std::vector<uint64_t> wigDist;
  std::vector<Peak> vPeak;
  double pthre;

 public:
  int32_t nbin;
  double nb_p, nb_n;
  WigStats(const int32_t _nbin=0, const double thre=0):
    ave(0), var(0), nb_p0(0),
    wigDist(WIGDISTNUM, 0),
    pthre(-log10(thre)),
    nbin(_nbin), nb_p(0), nb_n(0)
  {}
  
  void setWigStats(const WigArray &wigarray);
  
  /*  void setZINBParam(const std::vector<int32_t> &ar) {
    MyStatistics::moment<int32_t> mm(ar, 0);
    ave = mm.getmean();
    var = mm.getvar();
    nb_p = var ? ave/var : 0;
    if (nb_p>=1) nb_p = 0.9;
    if (nb_p<=0) nb_p = 0.1; 
    nb_n = ave * nb_p /(1 - nb_p);
    //    std::cout << ave << "\t" << var << "\t" << nb_p << "\t" << nb_n << std::endl;
    //    if(ave) estimateZINB(nb_p, nb_n);
    }*/

  /*  void estimateZINB(const double nb_p_pre, const double nb_n_pre) {
    uint32_t thre = getWigDistThre(wigDist, nbin);
    double parray[thre+1];
    parray[0] = thre;
    for(uint32_t i=0; i<thre; ++i) parray[i+1] = getpWig(i);
    iterateZINB(&parray, nb_p_pre, nb_n_pre, nb_p, nb_n, nb_p0);
    return;
    }*/

  void addWigDist(const WigStats &chr) {
    for (size_t i=0; i<wigDist.size(); ++i) wigDist[i] += chr.wigDist[i];
  }
    
  double getpWig(const int32_t i) const { return getratio(wigDist[i], nbin); }

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
  void printWigDist(std::ofstream &out, const int32_t i) const {
    out << boost::format("%1%\t%2%\t") % wigDist[i] % getpWig(i);
  }
  void printPoispar(std::ofstream &out) const {
    out << boost::format("%1$.3f\t%2$.3f\t") % ave % var;
  }
  void printZINBpar(std::ofstream &out) const {
    out << boost::format("%1%\t%2%\t%3%") % nb_p % nb_n % nb_p0;
  }
  int32_t getnbin() const { return nbin; }
  int32_t getWigDistsize() const { return wigDist.size(); }

  void peakcall(const WigArray &wigarray, const std::string chrname);

  int32_t printPeak(std::ofstream &out, const int32_t num, const int32_t binsize) const {
    for(uint32_t i=0; i<vPeak.size(); ++i) {
      vPeak[i].print(out, num+i, binsize);
    }
    return vPeak.size();
  }
  double getpthre() const { return pthre; }
  double getlogp(const double val) const { return getlogpZINB(val, nb_p, nb_n); }
};

class WigArray {
  std::vector<int32_t> array;
  double geta;
  
  template <class T> double rmGeta(const T val)  const { return val/geta; } 
  template <class T> double addGeta(const T val) const { return val*geta; }
  
  void checki(const size_t i) const {
    if(i>=array.size()) PRINTERR("Invalid i for WigArray: " << i << " > " << array.size());
  }

 public:
  WigArray(): geta(1000.0){}
  WigArray(const size_t num, const int32_t val):
    array(num, val), geta(1000.0)
  {}

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

  void Smoothing(const int32_t nsmooth) {
    GaussianSmoothing(array, nsmooth);
  }

  int64_t getArraySum() const {
    int64_t sum(0);
    for(auto x: array) sum += x;
    return rmGeta(sum);
  }
  double getMinValue() const {
    int32_t min(*std::min_element(array.begin(), array.end()));
    return rmGeta(min);
  }
  double getPercentile(double per) const {
    int32_t v95(MyStatistics::getPercentile(array, per));
    return rmGeta(v95);
  }

  void outputAsWig(FILE *File, const int32_t binsize, const int32_t showzero) const {
    for(size_t i=0; i<array.size(); ++i) {
      if(array[i] || showzero) fprintf(File, "%zu\t%.0f\n", i*binsize +1, rmGeta(array[i]));
    }
  }
  void outputAsBedGraph(FILE *File, const int32_t binsize, const std::string &name, const uint64_t chrend, const int32_t showzero) {
    for(size_t i=0; i<array.size()-1; ++i) {
      if(array[i] || showzero) fprintf(File, "%s %zu %zu %.0f\n", name.c_str(), i*binsize, (i+1) * binsize, rmGeta(array[i]));
    }
    size_t i = array.size()-1;
    if(array[i] || showzero) fprintf(File, "%s %zu %zu %.0f\n", name.c_str(), i*binsize, chrend, rmGeta(array[i]));
    
  }
  /*  void outputAsBinary(std::ofstream &out) const {
    for(size_t i=0; i<array.size(); ++i) out.write((char *)&array[i], sizeof(int32_t));
  }
  void readBinary(std::ifstream &in, const int32_t nbin) const {
    for(int32_t i=0; i<nbin; ++i) in.read((char *)&array[i], sizeof(int32_t));
    }*/
  void dump() const {
    for (auto x: array) std::cout << x << std::endl;
  }
};

class WigStatsGenome {
  int32_t binsize;
  int32_t rcenter;
  WigType type;
  int32_t outputzero;
  double pthre_inter;

public:
  std::vector<WigStats> chr;
  WigStats genome;

  WigStatsGenome(){}

  void setOpts(MyOpt::Opts &allopts) {
    MyOpt::Opts opt("Wigarray", 100);
    opt.add_options()
      ("outputformat",
       boost::program_options::value<int32_t>()->default_value(3)->notifier(boost::bind(&MyOpt::range<int32_t>, _1, 0, static_cast<int>(WigType::WIGTYPENUM) -2, "--outputformat")),
       "Output format\n   0: compressed wig (.wig.gz)\n   1: uncompressed wig (.wig)\n   2: bedGraph (.bedGraph)\n   3: bigWig (.bw)")
      ("binsize,b",
       boost::program_options::value<int32_t>()->default_value(100)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--binsize")),
       "bin size")
      ("outputzero",
       boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 0, "--outputzero")),
       "output zero-value bins (default: omitted)")
      ("rcenter",
       boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 0, "--rcenter")),
       "consider length around the center of fragment")
      ("pthre", boost::program_options::value<double>()->default_value(1e-3), "p-value threshold for peak calling")
      ;
    allopts.add(opt);
  }
  void setValues(const MyOpt::Variables &values, const std::vector<SeqStats> &_chr) {
    binsize = MyOpt::getVal<int32_t>(values, "binsize");
    rcenter = MyOpt::getVal<int32_t>(values, "rcenter");
    type    = static_cast<WigType>(MyOpt::getVal<int32_t>(values, "outputformat"));
    outputzero  = MyOpt::getVal<int32_t>(values, "outputzero");
    pthre_inter = MyOpt::getVal<double>(values, "pthre");

    for (auto &x: _chr) chr.emplace_back(x.getlen()/binsize +1, pthre_inter);
    for (auto &x: chr) genome.nbin += x.getnbin();
  }
  void dump() const {
    std::vector<std::string> strType = {"COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};
    std::cout << "\tOutput format: " << strType[static_cast<int32_t>(type)] << std::endl;
    std::cout << "Binsize: " << binsize << " bp" << std::endl;
    std::cout << "Peak calling threshold: " << pthre_inter << std::endl;
  }

  int32_t getbinsize() const { return binsize; }
  int32_t getWigDistsize() const { return genome.getWigDistsize(); }
  int32_t getrcenter() const { return rcenter; }
  int32_t isoutputzero() const { return outputzero; }
  WigType getWigType() const { return type; }

  void setWigStats(const int32_t id, const WigArray &array) {
    chr[id].setWigStats(array);
    genome.addWigDist(chr[id]);
  }
  /*  void estimateZINB(const int32_t id) {
    genome.estimateZINB(chr[id].nb_p, chr[id].nb_n);
   }*/
  void printPeak(const std::string &prefix) const {
    std::string filename = prefix + ".peak.xls";
    std::ofstream out(filename);

    Peak v;
    v.printHead(out);
    int32_t num(0);
    for (auto &x: chr) num += x.printPeak(out, num, binsize);
  }
};

#endif /* _WIGSTATS_HPP_ */
