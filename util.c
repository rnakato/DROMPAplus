/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "util.h"

void addmp(std::map<int, double> &mpto, const std::map<int, double> &mpfrom, double w)
{
  for(auto itr = mpfrom.begin(); itr != mpfrom.end(); ++itr) {
    mpto[itr->first] += itr->second * w;
  }
}

int chrname2int(std::string str)
{
  std::string chr;
  int chrnum(0);
  if(!str.find("chr")) chr = str.substr(3);
  else chr = str;
  try {
    chrnum = stoi(chr);
  } catch (std::invalid_argument e) {  // 数値以外
    chrnum = 0;
  }
  return chrnum;
}
