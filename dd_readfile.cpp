/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */

#include <unordered_map>
#include "dd_readfile.h"

/* 1:ChIP   2:Input   3:name   4:peaklist   5:binsize
   6:scale_tag   7:scale_ratio   8:scale_pvalue */
void scan_samplestr(string str, unordered_map<string, SampleFile> &sample, vector<SamplePair> &samplepair)
{
vector<string> v;
  boost::split(v, str, boost::algorithm::is_any_of(","));
  
  if(v.size() >8) {
    cerr << "error: sample string has ',' more than 8: " << str << endl;
    exit(1);
  }
  printf("tes3\n");
  
  if(v[0] == "") {
      cerr << "please specify ChIP sample: " << str << endl;
      exit(1);
  } else {
    if(sample.find(v[0]) == sample.end()) sample[v[0]] = SampleFile(v[0]);
  }
  if(v[1] != "") {
    if(sample.find(v[1]) == sample.end()) sample[v[1]] = SampleFile(v[1]);
  }

SamplePair p(v);
#ifdef DEBUG
  p.print();
#endif
  samplepair.push_back(p);
  
  return;
}
