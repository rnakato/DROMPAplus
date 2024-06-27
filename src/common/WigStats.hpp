/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _WIGSTATS_HPP_
#define _WIGSTATS_HPP_

#include <vector>
#include <fstream>
#include <boost/bind.hpp>
#include "extendBedFormat.hpp"
#include "statistics.hpp"
#include "util.hpp"
#include "../submodules/SSP/common/BoostOptions.hpp"
#include "SeqStatsDROMPA.hpp"

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

class WigArray {
  std::vector<int64_t> array;
  double geta;
  int32_t length_bin;
  enum {BINNUM_FOR_LOCALPOISSON=1000};

  template <class T> double rmGeta(const T val)  const { return val/geta; }
  template <class T> double addGeta(const T val) const { return val*geta; }

  void checki(const size_t i) const {
    if (i>=array.size())
      PRINTERR_AND_EXIT("Invalid i for WigArray: " << i << " > " << array.size());
  }

 public:
  WigArray(): geta(10000.0), length_bin(BINNUM_FOR_LOCALPOISSON) {}
  WigArray(const size_t num, const int32_t val, const int32_t _length_bin):
    array(num, val), geta(10000.0), length_bin(_length_bin)
  {
    if (length_bin <= 0) length_bin = BINNUM_FOR_LOCALPOISSON;
  }

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
  void divideval(const size_t i, const double val) {
    checki(i);
    array[i] /= val;
  }

  void Smoothing(const int32_t nsmooth, const int32_t smoothingtype) {
    if (smoothingtype == 1) GaussianSmoothing(array, nsmooth);
    else AverageSmoothing(array, nsmooth);
  }

  int64_t getArraySum() const {
    int64_t sum(0);
    for (auto x: array) sum += x;
    return rmGeta(sum);
  }

  int32_t getLengthBin() const { return length_bin; }

  double getMinValue() const {
    int32_t min(*std::min_element(array.begin(), array.end()));
    return rmGeta(min);
  }

  double getLocalAverage(const int32_t i) const {
//    int32_t length_bin(BINNUM_FOR_LOCALPOISSON);
    checki(i);
    if (i<0) PRINTERR_AND_EXIT("Invalid i for WigArray: " << i << " < 0.");

   // printf("length_bin = %d\n", length_bin);

    int64_t ave(0);
    int32_t lenhalf(length_bin/2);
    int32_t left(std::max(i-lenhalf, 0));
    int32_t right(std::min(i+lenhalf, (int32_t)array.size()));
    for (int32_t j=left; j<right; ++j) ave += array[j];
    ave /= right - left;
    return rmGeta(ave);
  }

  double getPercentile(double per) const {
    int32_t v95(MyStatistics::getPercentile(array, per));
    return rmGeta(v95);
  }

  void outputAsWig(FILE *File, const int32_t binsize, const int32_t showzero, const bool isfloat) const {
    for (size_t i=0; i<array.size(); ++i) {
      if (array[i] || showzero) {
	if (isfloat) fprintf(File, "%zu\t%.3f\n", i*binsize +1, rmGeta(array[i]));
	else         fprintf(File, "%zu\t%.0f\n", i*binsize +1, rmGeta(array[i]));
      }
    }
  }
  void outputAsBedGraph(FILE *File, const int32_t binsize, const std::string &name, const uint64_t chrend, const int32_t showzero, const bool isfloat) {
    for (size_t i=0; i<array.size()-1; ++i) {
      if (array[i] || showzero) {
	if (isfloat) fprintf(File, "%s\t%zu\t%zu\t%.3f\n", name.c_str(), i*binsize, (i+1) * binsize, rmGeta(array[i]));
	else         fprintf(File, "%s\t%zu\t%zu\t%.0f\n", name.c_str(), i*binsize, (i+1) * binsize, rmGeta(array[i]));
      }
    }
    size_t i = array.size()-1;
    if (array[i] || showzero) {
      if (isfloat) fprintf(File, "%s\t%zu\t%lu\t%.3f\n", name.c_str(), i*binsize, (uint64_t)chrend, rmGeta(array[i]));
      else         fprintf(File, "%s\t%zu\t%lu\t%.0f\n", name.c_str(), i*binsize, (uint64_t)chrend, rmGeta(array[i]));
    }
  }
  void dump() const {
    for (auto &x: array) std::cout << x << std::endl;
  }

};

class WigStats {
  enum {WIGDISTNUM=200};

  std::vector<uint64_t> wigDist;

 public:
  int32_t nbin;
  int32_t binsize;
  double nb_p, nb_n;
  WigStats(const int32_t _nbin=0, const int32_t _binsize=0):
    wigDist(WIGDISTNUM, 0),
    nbin(_nbin), binsize(_binsize),
    nb_p(0), nb_n(0)
  {}

  void setWigStats(const WigArray &wigarray);

  void addWigDist(const WigStats &chr) {
    for (size_t i=0; i<wigDist.size(); ++i) wigDist[i] += chr.wigDist[i];
  }
  void printWigDist(std::ofstream &out, const int32_t i) const {
    out << boost::format("%1%\t%2%\t") % wigDist[i] % getpWig(i);
  }
  double getpWig(const int32_t i) const { return getratio(wigDist[i], nbin); }

/*  double getPoisson(const int32_t i) const {
    if (ave) return _getPoisson(i, ave);
    else return 0;
  }*/
/*  double getNegativeBinomial(const int32_t i) const {
    return _getNegativeBinomial(i, nb_p, nb_n);
  }
  double getZINB(const int32_t i) const {
    if (ave) return _getZINB(i, nb_p, nb_n, nb_p0);
    else return 0;
  }
  void printPoispar(std::ofstream &out) const {
    out << boost::format("%1$.3f\t%2$.3f\t") % ave % var;
  }
  void printZINBpar(std::ofstream &out) const {
    out << boost::format("%1%\t%2%\t%3%") % nb_p % nb_n % nb_p0;
  }*/
  int32_t getnbin() const { return nbin; }
  int32_t getWigDistsize() const { return wigDist.size(); }

  /*  void setZINBParam(const std::vector<int32_t> &ar) {
    MyStatistics::moment<int32_t> mm(ar, 0);
    ave = mm.getmean();
    var = mm.getvar();
    nb_p = var ? ave/var : 0;
    if (nb_p>=1) nb_p = 0.9;
    if (nb_p<=0) nb_p = 0.1;
    nb_n = ave * nb_p /(1 - nb_p);
    //    std::cout << ave << "\t" << var << "\t" << nb_p << "\t" << nb_n << std::endl;
    //    if (ave) estimateZINB(nb_p, nb_n);
    }*/

  /*  void estimateZINB(const double nb_p_pre, const double nb_n_pre) {
    uint32_t thre = getWigDistThre(wigDist, nbin);
    double parray[thre+1];
    parray[0] = thre;
    for (uint32_t i=0; i<thre; ++i) parray[i+1] = getpWig(i);
    iterateZINB(&parray, nb_p_pre, nb_n_pre, nb_p, nb_n, nb_p0);
    return;
    }*/
};

class WigStatsGenome {
  int32_t binsize;
  int32_t rcenter;
  WigType type;
  bool outputzero;
  bool onlyreadregion;

public:
  std::vector<WigStats> chr;
  WigStats genome;

  WigStatsGenome(): binsize(0), rcenter(0), type(WigType::NONE), outputzero(false), onlyreadregion(false) {}

  void setOpts(MyOpt::Opts &allopts) {
    MyOpt::Opts opt("Wigarray", 100);
    opt.add_options()
      ("outputformat",
       boost::program_options::value<int32_t>()->default_value(3)->notifier(boost::bind(&MyOpt::range<int32_t>, _1, 0, static_cast<int>(WigType::WIGTYPENUM) -2, "--outputformat")),
       "Output format\n   0: compressed wig (.wig.gz)\n   1: uncompressed wig (.wig)\n   2: bedGraph (.bedGraph)\n   3: bigWig (.bw)")
      ("binsize,b",
       boost::program_options::value<int32_t>()->default_value(100)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--binsize")),
       "bin size")
      ("outputzero", "output zero-value bins (default: omitted)")
      ("rcenter",
       boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 0, "--rcenter")),
       "consider length around the center of fragment")
      ("onlyreadregion", "(for paired-end) count only read region (default: full fragment length)")
      ;
    allopts.add(opt);
  }
  void setValues(const MyOpt::Variables &values, const std::vector<SeqStats> &_chr) {
    binsize = MyOpt::getVal<int32_t>(values, "binsize");
    rcenter = MyOpt::getVal<int32_t>(values, "rcenter");
    type    = static_cast<WigType>(MyOpt::getVal<int32_t>(values, "outputformat"));
    outputzero = values.count("outputzero");
    onlyreadregion = values.count("onlyreadregion");

    for (auto &x: _chr) chr.emplace_back(x.getlen()/binsize +1);
    for (auto &x: chr) genome.nbin += x.getnbin();
  }
  void dump() const {
    std::vector<std::string> strType = {"COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};
    std::cout << "Output format: " << strType[static_cast<int32_t>(type)] << std::endl;
    std::cout << "Binsize: " << binsize << " bp" << std::endl;
  }

  int32_t getbinsize() const { return binsize; }
  int32_t getWigDistsize() const { return genome.getWigDistsize(); }
  int32_t getrcenter() const { return rcenter; }
  bool isoutputzero() const { return outputzero; }
  bool isonlyreadregion() const { return onlyreadregion; }
  WigType getWigType() const { return type; }

  void setWigStats(const int32_t id, const WigArray &array) {
    chr[id].setWigStats(array);
    genome.addWigDist(chr[id]);
  }
  /*  void estimateZINB(const int32_t id) {
    genome.estimateZINB(chr[id].nb_p, chr[id].nb_n);
   }*/
};

#endif /* _WIGSTATS_HPP_ */
