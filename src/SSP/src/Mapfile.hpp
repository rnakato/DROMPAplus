/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _MAPFILE_HPP_
#define _MAPFILE_HPP_

#include <boost/format.hpp>
#include "SeqStats.hpp"
#include "MThread.hpp"
#include "../common/BedFormat.hpp"
#include "../common/util.hpp"
#include "../common/BoostOptions.hpp"

class FragmentLengthDist {
  MyOpt::Opts opt;

  enum {ReadLenMax=200, FragLenMax=1000};
  int32_t flen_ssp;
  int32_t flen_def;
  std::vector<int32_t> vlenF3;
  std::vector<int32_t> vlenF5;
  std::vector<int32_t> vflen;
  int32_t nomodel;
  int32_t pairedend;

  template <class T>
  void printVector(std::ofstream &out, const std::vector<T> v, const std::string &str, const uint64_t nread)
  {
    out << str << " length distribution" << std::endl;
    out << "length\tnumber\tproportion" << std::endl;
    for(size_t i=0; i<v.size(); ++i)
      if(v[i]) out << boost::format("%1%\t%2%\t%3%\n") % i % v[i] % getratio(v[i], nread);
  }
 
 public:
 FragmentLengthDist():
   opt("Fragment",100),
   flen_ssp(0),
   vlenF3(ReadLenMax,0),
   vlenF5(ReadLenMax,0),
   vflen(FragLenMax,0)
  {
    opt.add_options()
      ("nomodel", "omit fraglent length estimation (default: estimated by strand-shift profile)")
      ("flen",
       boost::program_options::value<int32_t>()->default_value(150)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--binsize")),
       "predefined fragment length (with --nomodel option)")
      ;
  }

  void setOpts(MyOpt::Opts &allopts) {
    allopts.add(opt);
  }
  
  void setValues(const MyOpt::Variables &values) {
    DEBUGprint("FragmentLengthDist setValues...");
    
    flen_def  = MyOpt::getVal<int32_t>(values, "flen");
    nomodel   = values.count("nomodel");
    pairedend = values.count("pair");
    
    DEBUGprint("FragmentLengthDist setValues done.");
  }

  int32_t isnomodel() const { return nomodel; }
  
  int32_t getlenF3 () const { return getmaxi(vlenF3); }
  int32_t getlenF5 () const { return getmaxi(vlenF5); }
  int32_t getflen4paired () const { return getmaxi(vflen); }

  void setflen_ssp(const int32_t len) { flen_ssp = len; }
  int32_t getflen() const {
    if(pairedend)     return getflen4paired();
    else if(!nomodel) return flen_ssp;
    else              return flen_def;
  }
  void printFlen(std::ofstream &out) const {
    if(pairedend)     out << "Most likely fragment length: " << getflen4paired() << std::endl;
    else if(!nomodel) out << "Estimated fragment length: "   << flen_ssp << std::endl;
    else              out << "Predefined fragment length: "  << flen_def << std::endl;
  }
  void addF3(const int32_t lenF3) { ++vlenF3[lenF3]; }
  void addF5(const int32_t lenF5) { ++vlenF5[lenF5]; }
  void addvflen(const int32_t flen) { ++vflen[flen]; }

  void outputDistFile(const std::string &prefix, const uint64_t nread);
};

class SeqStatsGenome {
  MyOpt::Opts opt;

  std::string inputfilename;
  std::string genometable;

  int32_t pairedend;
  int32_t maxins;
  int32_t specifyFtype;
  std::string ftype;

  std::string name;
  double depth;
  double sizefactor;

  void readGenomeTable(const std::string &gt);

 public:
  std::vector<SeqStats> chr;
  std::vector<MyMthread::chrrange> vsepchr;
  FragmentLengthDist dflen;
  
 SeqStatsGenome(): 
   opt("Genome",100), ftype(""),
   name("Genome"), depth(0), sizefactor(0) {
   using namespace boost::program_options;
   opt.add_options()
     ("gt", value<std::string>(),
      "Genome table (tab-delimited file describing the name and length of each chromosome)")
     ("mptable", value<std::string>(),
      "Genome table of mappable regions")
      ;
 }

  const std::string & getInputfile() const { return inputfilename; }
  const std::string & getGenomeTable() const { return genometable; }
  int32_t isPaired()    const { return pairedend; }
  int32_t getmaxins()   const { return maxins; }
  int32_t onFtype()     const { return specifyFtype; }
  const std::string & getftype() const { return ftype; }

  void setOpts(MyOpt::Opts &allopts) {
    allopts.add(opt);
    dflen.setOpts(allopts);
  }
  void setValues(const MyOpt::Variables &values);

  std::string getname() const { return name; }
  uint64_t getlen() const {
    uint64_t len(0);
    for(auto &x:chr) len += x.getlen();
    return len;
  }
  uint64_t getlenmpbl() const {
    uint64_t len_mpbl(0);
    for(auto &x:chr) len_mpbl += x.getlenmpbl();
    return len_mpbl;
  }
  double getpmpbl() const {
    return getratio(getlenmpbl(), getlen());
  }
  uint64_t getnread (const Strand::Strand strand) const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread(strand);
    return nread;
  }
  uint64_t getnread_nonred (const Strand::Strand strand) const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread_nonred(strand);
    return nread;
  }
  uint64_t getnread_red (const Strand::Strand strand) const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread_red(strand);
    return nread;
  }
  uint64_t getnread_rpm (const Strand::Strand strand) const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread_rpm(strand);
    return nread;
  }
  uint64_t getnread_afterGC (const Strand::Strand strand) const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread_afterGC(strand);
    return nread;
  }
  uint64_t getnread_inbed() const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread_inbed();
    return nread;
  }
  double getFRiP() const {
    return getratio(getnread_inbed(), getnread_nonred(Strand::BOTH));
  }
  void setdepth(const double d) { depth = d; }
  double getdepth() const { return depth; }
  double getsizefactor()const { return sizefactor; }

  void setsizefactor(const double w) { sizefactor = w; }

  void printReadstats() const {
    std::cout << "name\tlength\tlen_mpbl\tread num\tnonred num\tred num\tnormed\tafterGC\tdepth" << std::endl;
    printSeqStats(*this);
    for(auto &x: chr) printSeqStats(x);
  }
};

int32_t setIdLongestChr(SeqStatsGenome &genome);

#endif /* _MAPFILE_HPP_ */
