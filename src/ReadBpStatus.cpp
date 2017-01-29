/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <iostream>
#include <string>
#include <cstdint>
#include <boost/algorithm/string.hpp>
#include "ReadBpStatus.hpp"
#include "SSP/common/inline.hpp"
#include "SSP/common/util.hpp"

std::vector<int32_t> readMpbl(std::string mpfile, std::string chrname, int32_t binsize, int32_t nbin)
{
  std::string filename = mpfile + "/map_fragL150_" + chrname + "_bin" + IntToString(binsize) +".txt";
  std::vector<int32_t> mparray(nbin, 0);

  isFile(filename);
  std::ifstream in(filename);

  std::string lineStr;
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty()) continue;
    std::vector<std::string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));

    int32_t n(stoi(v[0])/binsize);
    double val(stof(v[1])*binsize);
    mparray[n] = val;
  }

  return mparray;
}

std::vector<BpStatus> readMpbl_binary(std::string mpfile, std::string chrname, int32_t chrlen)
{
  std::string filename = mpfile + "/map_" + chrname + "_binary.txt";
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

std::vector<BpStatus> readMpbl_binary(int32_t chrlen)
{
  std::vector<BpStatus> mparray(chrlen, BpStatus::MAPPABLE);
  return mparray;
}

std::vector<BpStatus> arraySetBed(std::vector<BpStatus> &array, std::string chrname, const std::vector<bed> &vbed)
{
  for(auto &bed: vbed) {
    if(bed.chr == chrname) {
      size_t s(bed.start);
      size_t e(bed.end);
      if(e >= array.size()) {
	std::cerr << "Warning: bedfile" << bed.start <<"-"<<bed.end << " > array size " << array.size() << std::endl;
	e = array.size()-1;
      }
      for(size_t i=s; i<=e; ++i) array[i] = BpStatus::INBED;
    } 
  }

  return array;
}
