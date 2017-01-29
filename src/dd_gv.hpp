/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_GV_H_
#define _DD_GV_H_

#include <iostream>
#include <unordered_map>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include "SSP/common/inline.hpp"
#include "SSP/common/seq.hpp"
#include "MyType.hpp"

class SampleFile {
  double lambda;
  double nb_p, nb_n, nb_p0;
  WigType iftype;
  int binsize;
 public:
    // *genome, *chr;
  std::vector<int> data;
 SampleFile() {}
 SampleFile(std::string &str) {
   std::vector<std::string> v;
   boost::split(v, str, boost::algorithm::is_any_of("."));
   int last(v.size()-1);
   if(v[last] == "wig") iftype = WigType::UNCOMPRESSWIG;    
   else if(v[last] == "gz" && v[last-1] == "wig") {
     iftype = WigType::COMPRESSWIG;
     --last;
   } else if(v[last] == "bedGraph") iftype = WigType::BEDGRAPH;
   else if(v[last] == "bin")        iftype = WigType::BINARY;
   else PRINTERR("invalid postfix: " << str);
   try {
     binsize = stoi(v[last-1]);
   }catch (...) {
     binsize = 0;
   }
   if(binsize <= 0) PRINTERR("invalid binsize: " << str);
   //   cout << iftype << "\t" << binsize << endl;
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
  std::string argvChIP, argvInput;
  //  double fc; //comp

  std::string peak_argv;
  /*  Peak *peak;
  char *peak_argv;
  char *peakarray;
  char *linename;
  int *binnum;*/
  
 public:
  std::string name;
  yScale scale;

 SamplePair(std::vector<std::string> v): argvInput(""), peak_argv(""), name("") {
    if(v[0] != "") argvChIP  = v[0];
    if(v.size() >=2 && v[1] != "") argvInput = v[1];
    if(v.size() >=3 && v[2] != "") name      = v[2];
    if(v.size() >=4 && v[3] != "") peak_argv = v[3];
    //    if(v.size() >=5 && v[4] != "") binsize   = stoi(v[4]);
    if(v.size() >=6 && v[5] != "") scale.tag = stof(v[5]);
    if(v.size() >=7 && v[6] != "") scale.ratio  = stof(v[6]);
    if(v.size() >=8 && v[7] != "") scale.pvalue = stof(v[7]);
  }
  void print() {
    std::cout << boost::format("ChIP: %1% label: %2% peaklist: %3%\n") % argvChIP % name % peak_argv;
    std::cout << boost::format("   Input: %1%\n") % argvInput;
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

class Param {
 public:
  std::vector<chrsize> gt;
  std::unordered_map<std::string, SampleFile> sample;
  std::vector<SamplePair> samplepair;
  std::vector<pdSample> pd;
  Param(){}
};

#endif /* _DD_GV_H_ */
