/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_GV_H_
#define _DD_GV_H_

#include <unordered_map>
#include "readdata.h"
#include "macro.h"

using namespace std;

class SampleFile {
  double lambda;
  double nb_p, nb_n, nb_p0;
 public:
  string argv;
    // *genome, *chr;
  vector<int> data;
 SampleFile(){}
 SampleFile(string str): argv(str) {}
};

class yScale {
 public:
  double tag;
  double ratio;
  double pvalue;
 yScale(): tag(0), ratio(0), pvalue(0) {}
};

class SamplePair {
  int overlay;
  string argvChIP, argvInput;
  //  double fc; //comp
  int binsize;

  string peak_argv;
  /*  Peak *peak;
  char *peak_argv;
  char *peakarray;
  char *linename;
  int *binnum;*/
  
 public:
  string name;
  yScale scale;

 SamplePair(vector<string> v): argvInput(""), binsize(0), peak_argv(""), name("") {
    if(v[0] != "") argvChIP  = v[0];
    if(v.size() >=2 && v[1] != "") argvInput = v[1];
    if(v.size() >=3 && v[2] != "") name      = v[2];
    if(v.size() >=4 && v[3] != "") peak_argv = v[3];
    if(v.size() >=5 && v[4] != "") binsize   = stoi(v[4]);
    if(v.size() >=6 && v[5] != "") scale.tag = stof(v[5]);
    if(v.size() >=7 && v[6] != "") scale.ratio  = stof(v[6]);
    if(v.size() >=8 && v[7] != "") scale.pvalue = stof(v[7]);
  }
  void print() {
    BPRINT("ChIP:%1% Input:%2% name:%3% peak:%4% binsize:%5% scale_tag %6% scale_ratio %7% scale_pvalue %8%\n")
      % argvChIP % argvInput % name % peak_argv % binsize % scale.tag % scale.ratio % scale.pvalue;
  }
};


#endif /* _DD_GV_H_ */
