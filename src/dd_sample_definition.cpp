/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "dd_sample_definition.hpp"

namespace {

}

SampleInfo::SampleInfo(const std::string &filename,
		       const std::vector<chrsize> &gt,
		       const int32_t b,
		       const WigType &type):
  binsize(0), totalreadnum(0), prefix("")
{
  std::vector<std::string> v;
  ParseLine(v, filename, '.');
  int32_t last(v.size()-1);

  if (type != WigType::NONE) iftype = type;
  else {
    if (v[last] == "wig") iftype = WigType::UNCOMPRESSWIG;
    else if (v[last] == "gz" && v[last-1] == "wig") {
      iftype = WigType::COMPRESSWIG;
      --last;
    } else if (v[last] == "bedGraph") iftype = WigType::BEDGRAPH;
    else if (v[last] == "bw")         iftype = WigType::BIGWIG;
    //     else if (v[last] == "bin")        iftype = WigType::BINARY;
    else PRINTERR_AND_EXIT("invalid postfix: " << filename);
  }
  setbinsize(v[last-1], b);
  for (int32_t i=0; i<last; ++i) prefix += v[i] + ".";
  gettotalreadnum(filename, gt);
}

void SampleInfo::setbinsize(std::string &v, const int32_t b)
{
  if (b>0) binsize = b;
  else {
    try {
      binsize = stoi(v);
    } catch (...) {
      binsize = 0;
    }
  }
  if (binsize <= 0) PRINTERR_AND_EXIT("invalid binsize: " << v);
}


SamplePairEach::SamplePairEach(const std::string &str, const vSampleInfo &vsinfo):
  binsize(0), argvChIP(""), argvInput(""), peak_argv(""), label(""), ratio(1)
{
  std::vector<std::string> v;
  ParseLine(v, str, ',');

  /* 1:ChIP   2:Input   3:label   4:peaklist   5:binsize
     6:scale_tag   7:scale_ratio   8:scale_pvalue */
  if (v[0] != "") argvChIP = v[0];
  if (v.size() >=2 && v[1] != "") argvInput = v[1];
  if (v.size() >=3 && v[2] != "") label     = v[2];
  if (v.size() >=4 && v[3] != "") peak_argv = v[3];
  if (peak_argv != "") peaks = parseBed_Hash<bed>(peak_argv);
  binsize = vsinfo.getbinsize(argvChIP);
  if (v.size() >=6 && v[5] != "") scale.tag = stod(v[5]);
  if (v.size() >=7 && v[6] != "") scale.ratio = stod(v[6]);
  if (v.size() >=8 && v[7] != "") scale.pvalue = stod(v[7]);

  //    printBed_Hash(peaks);
}

std::vector<bed> SamplePairEach::getpeaksChr(const std::string &chrname) const
{
  if (peak_argv != "") return peaks.at(rmchr(chrname));
  else return std::vector<bed>();
}

void SamplePairEach::print() const
{
  std::cout << boost::format("ChIP: %1% label: %2% peaklist: %3%\n") % argvChIP % label % peak_argv;
  std::cout << boost::format("   Input: %1%\n") % argvInput;
  std::cout << boost::format("   binsize: %1%\n") % binsize;
}
