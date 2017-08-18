/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <unordered_map>
#include "dd_readfile.hpp"
#include "dd_gv.hpp"

/* 1:ChIP   2:Input   3:name   4:peaklist   5:binsize
   6:scale_tag   7:scale_ratio   8:scale_pvalue */
void scan_samplestr(const std::string &str, std::unordered_map<std::string, SampleFile> &sample, std::vector<SamplePair> &samplepair)
{
std::vector<std::string> v;
  boost::split(v, str, boost::algorithm::is_any_of(","));

  if(v.size() >8) {
    std::cerr << "error: sample std::string has ',' more than 8: " << str << std::endl;
    exit(1);
  }
  
  if(v[0] == "") {
      std::cerr << "please specify ChIP sample: " << str << std::endl;
      exit(1);
  } else {
    if(sample.find(v[0]) == sample.end()) sample[v[0]] = SampleFile(v[0]);
  }
  if(v.size() >=2 && v[1] != "") {
    if(sample.find(v[1]) == sample.end()) sample[v[1]] = SampleFile(v[1]);
    if(sample[v[0]].getbinsize() != sample[v[1]].getbinsize()) PRINTERR("binsize of ChIP and Input should be same. " << str);
  }

  samplepair.emplace_back(v);
  
  return;
}

pdSample scan_pdstr(const std::string &str)
{
  std::vector<std::string> v;
  boost::split(v, str, boost::algorithm::is_any_of(","));

  if(v.size() >2) {
    std::cerr << "error: sample std::string has ',' more than 2: " << str << std::endl;
    exit(1);
  }

  pdSample pd;
  if(v[0] == "") {
      std::cerr << "please specify file: " << str << std::endl;
      exit(1);
  } else {
    pd.argv = v[0];
  }
  if(v[1] != "") pd.name = v[1];
  else pd.name = v[0];

  return pd;
}
