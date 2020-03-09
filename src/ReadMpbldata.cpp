/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <iostream>
#include <fstream>
#include <vector>
#include "ReadMpbldata.hpp"

std::vector<int32_t> readMpbl(const std::string &mpfile,
			      const std::string &chrname,
			      const int32_t binsize,
			      const int32_t nbin)
{
  std::string filename = mpfile + "/map_fragL150_" + chrname + "_bin" + std::to_string(binsize) +".txt";
  std::vector<int32_t> mparray(nbin, 0);

  isFile(filename);
  std::ifstream in(filename);

  std::string lineStr;
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty()) continue;
    std::vector<std::string> v;
    ParseLine(v, lineStr, '\t');

    int32_t n(stoi(v[0])/binsize);
    double val(stod(v[1])*binsize);
    mparray[n] = val;
  }

  return mparray;
}

std::vector<BpStatus> readMpbl_binary(const std::string &mpfile,
				      const std::string &chrname,
				      const int32_t chrlen)
{
  static int on(0);

  if(mpfile == "") {
    if(!on) {
      std::cout << "Mappability file is not specified. All genomeic regions are considered as mappable." << std::endl;
      on=1;
    }
    return std::vector<BpStatus>(chrlen, BpStatus::MAPPABLE);
  }

  std::string filename = mpfile + "/map_" + chrname + "_binary.txt";

  if(!on) {
    std::cout << "Mappability file is specified: " << filename << std::endl;
    on=1;
  }
  std::vector<BpStatus> mparray(chrlen, BpStatus::UNMAPPABLE);

  isFile(filename);
  int32_t n(0);
  int8_t c;
  std::ifstream in(filename);
  while (!in.eof()) {
    c = in.get();
    if(c==' ') continue;
    if(c=='1') mparray[n] = BpStatus::MAPPABLE;
    ++n;
    if(n >= chrlen-1) break;
  }

  return mparray;
}

void OverrideBedToArray(std::vector<BpStatus> &array,
			const std::string &chrname,
			const std::vector<bed> &vbed)
{
  int32_t chrlen(array.size());
  for(auto &bed: vbed) {
    if(bed.chr == chrname) {
      size_t s(std::max(0, bed.start));
      size_t e(std::min(bed.end, chrlen-1));
      for(size_t i=s; i<=e; ++i) array[i] = BpStatus::INBED;
    }
  }
  return;
}
