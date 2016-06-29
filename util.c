/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "util.h"

std::string rmchr(const std::string &chr)
{
  std::string s;
  if(!chr.find("chr")) s = chr.substr(3);
  else s = chr;
  return s;
}

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
    if(chr=="I")         chrnum = 1;
    else if(chr=="II")   chrnum = 2;
    else if(chr=="III")  chrnum = 3;
    else if(chr=="IV")   chrnum = 4;
    else if(chr=="V")    chrnum = 5;
    else if(chr=="VI")   chrnum = 6;
    else if(chr=="VII")  chrnum = 7;
    else if(chr=="VIII") chrnum = 8;
    else if(chr=="IX")   chrnum = 9;
    else if(chr=="X")    chrnum = 10;
    else if(chr=="XI")   chrnum = 11;
    else if(chr=="XII")  chrnum = 12;
    else if(chr=="XIII") chrnum = 13;
    else if(chr=="XIV")  chrnum = 14;
    else if(chr=="XV")   chrnum = 15;
    else if(chr=="XVI")  chrnum = 16;
    chrnum = 0;
  }
  return chrnum;
}
