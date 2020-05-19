/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _SEQSTATSDROMPA_HPP_
#define _SEQSTATSDROMPA_HPP_

#include "../submodules/SSP/src/SeqStats.hpp"

class SeqStats : public SeqStatsSSP {
 public:

 SeqStats(std::string &s, int32_t l):
   SeqStatsSSP(s, l) {}

  void setFRiP(const std::vector<bed> &vbed);
  
  double getFRiP() const {
    return getratio(nread_inbed, getnread_nonred(Strand::BOTH));
  }
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
   opt("Genome",100),
   pairedend(0), maxins(0), specifyFtype(0),
   ftype(""),
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


#endif /* _SEQSTATSDROMPA_HPP_ */
