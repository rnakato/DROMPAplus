/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_CLASS_H_
#define _DD_CLASS_H_

#include <boost/algorithm/string.hpp>
#include "WigStats.hpp"

class chrsize;

class SampleFile {
  double lambda;
  double nb_p, nb_n, nb_p0;
  WigType iftype;
  int32_t binsize;
 public:
  std::vector<int> data;

  void setbinsize(std::string &v, const int32_t b) {
    if(b>0) binsize = b;
    else {
      try {
	binsize = std::stoi(v);
      }catch (...) {
	binsize = 0;
      }
    }
    if(binsize <= 0) PRINTERR("invalid binsize: " << v);
  }
  
 SampleFile() {}
  SampleFile(const std::string &str, const int32_t b, const WigType &type): binsize(0) {
   std::vector<std::string> v;
   boost::split(v, str, boost::algorithm::is_any_of("."));
   int last(v.size()-1);

   // v.size()==1の場合分け

   if(type != WigType::NONE) {
     iftype = type;
     binsize = b;
   } else {
     if(v[last] == "wig") iftype = WigType::UNCOMPRESSWIG;    
     else if(v[last] == "gz" && v[last-1] == "wig") {
       iftype = WigType::COMPRESSWIG;
       --last;
     } else if(v[last] == "bedGraph") iftype = WigType::BEDGRAPH;
     else if(v[last] == "bw")         iftype = WigType::BIGWIG;
     else if(v[last] == "bin")        iftype = WigType::BINARY;
     else PRINTERR("invalid postfix: " << str);
   }
   setbinsize(v[last-1], b);
  }
 int getbinsize()    const { return binsize; }
 WigType getiftype() const { return iftype; }
};

class yScale {
  enum {TAG_DEFAULT=30, RATIO_DEFAULT=3, P_DEFAULT=5};
 public:
  double tag;
  double ratio;
  double pvalue;
 yScale(): tag(TAG_DEFAULT), ratio(RATIO_DEFAULT), pvalue(P_DEFAULT) {}
};

class SamplePair {
  int overlay;
  //  double fc; //comp

  std::string peak_argv;
  /*  Peak *peak;
  char *peak_argv;
  char *peakarray;
  char *linename;
  int *binnum;*/
  int binsize;
  
 public:
  std::string argvChIP, argvInput;
  std::string name;
  yScale scale;

  SamplePair(const std::vector<std::string> &v, const int32_t b): peak_argv(""), argvInput(""), name("") {
    if(v[0] != "") argvChIP  = v[0];
    if(v.size() >=2 && v[1] != "") argvInput = v[1];
    if(v.size() >=3 && v[2] != "") name      = v[2];
    if(v.size() >=4 && v[3] != "") peak_argv = v[3];
    binsize = b;
    if(v.size() >=6 && v[5] != "") scale.tag = stof(v[5]);
    if(v.size() >=7 && v[6] != "") scale.ratio  = stof(v[6]);
    if(v.size() >=8 && v[7] != "") scale.pvalue = stof(v[7]);
  }
  void print() {
    std::cout << boost::format("ChIP: %1% label: %2% peaklist: %3%\n") % argvChIP % name % peak_argv;
    std::cout << boost::format("   Input: %1%\n") % argvInput;
    std::cout << boost::format("   binsize: %1%\n") % binsize;
  }
  void printall() {
    std::cout << boost::format("ChIP:%1% Input:%2% name:%3% peak:%4% scale_tag %5% scale_ratio %6% scale_pvalue %7%\n")
      % argvChIP % argvInput % name % peak_argv % scale.tag % scale.ratio % scale.pvalue;
  }
};

class pdSample {
 public:
  std::string argv;
  std::string name;
  pdSample(){}
};


namespace DROMPA {
  class Global {
    WigType iftype;
    
    // smoothing
  public:
    bool includeYM;
    std::vector<chrsize> gt;
    std::unordered_map<std::string, SampleFile> sample;
    std::vector<SamplePair> samplepair;
    std::vector<pdSample> pd;
    
    Global(): iftype(WigType::NONE), includeYM(false) {}
    void getWigType(const int32_t t) {
      iftype = static_cast<WigType>(t);
    }

    WigType getIfType() const {return iftype;}
  };

}


#endif /* _DD_CLASS_H_ */
