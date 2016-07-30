/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef READGENE_H
#define READGENE_H

#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "seq.h"
#include "macro.h"

using GeneDataMap = std::unordered_map<std::string, genedata>;
using HashOfGeneDataMap = std::unordered_map<std::string, GeneDataMap>;

enum bpstatus {UNMAPPABLE, INBED, MAPPABLE, COVREAD_ALL, COVREAD_NORM};

int countmp(HashOfGeneDataMap &);
std::vector<std::string> scanGeneName(const HashOfGeneDataMap &);
HashOfGeneDataMap extract_mp(const HashOfGeneDataMap &, const std::vector<std::string>);
std::vector<std::string> readGeneList(const std::string&);
HashOfGeneDataMap parseRefFlat(const std::string&);
HashOfGeneDataMap parseGtf(const std::string&, const int);
HashOfGeneDataMap construct_gmp(const HashOfGeneDataMap &);
void printMap(const HashOfGeneDataMap &);
void printRefFlat(const HashOfGeneDataMap &);
std::vector<chrsize> read_genometable(const std::string&);

std::vector<int>  readMpbl(std::string, std::string, int, int);
std::vector<char> readMpbl_binary(int);
std::vector<char> readMpbl_binary(std::string, std::string, int);
std::vector<char> arraySetBed(std::vector<char> &, std::string, std::vector<bed>);
std::string IntToString(int n);

template <class T>
std::vector<T> parseBed(const std::string &fileName)
{
  std::vector<T> vbed;
  std::ifstream in(fileName);
  if(!in) PRINTERR("BED file does not exist.");

  std::string lineStr;
  std::vector<std::string> v;
  while (!in.eof()) {
    getline(in, lineStr);

    if(lineStr.empty() || lineStr[0] == '#' || !lineStr.find("chromosome")) continue;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    T bed(v);
    vbed.push_back(bed);
  }

  return vbed;
}

template <class T>
void printBed(const std::vector<T> &vbed)
{
  for (auto x: vbed) {
    x.print();
    std::cout << std::endl;
  }
  std::cout << "bed num: " << vbed.size() << std::endl;
  return;
}

#endif  // READGENE_H
