#ifndef READGENE_H
#define READGENE_H

#include <iostream>
#include <unordered_map>
#include <map>
#include <string>
#include <fstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include "seq.h"
#include "warn.h"

using namespace std;

int countmp(unordered_map<string, unordered_map<string, genedata>> &);
vector<string> scanGeneName(const unordered_map<string, unordered_map<string, genedata>> &);
unordered_map<string, unordered_map<string, genedata>> extract_mp(const unordered_map<string, unordered_map<string, genedata>> &, const vector<string>);
vector<string> readGeneList(const string&);
unordered_map<string, unordered_map<string, genedata>> parseRefFlat(const string&);
unordered_map<string, unordered_map<string, genedata>> parseGtf(const string&, const int);
unordered_map<string, unordered_map<string, genedata>> construct_gmp(const unordered_map<string, unordered_map<string, genedata>> &);
void printMap(const unordered_map<string, unordered_map<string, genedata>> &);
void printRefFlat(const unordered_map<string, unordered_map<string, genedata>> &);
map<string, int> read_genometable(const string&);
  
template <class T>
vector<T> parseBed(const string& fileName)
{
  vector<T> vbed;
  ifstream in(fileName);
  if(!in) printerr("BED file does not exist.");
  
  string lineStr;
  vector<string> v;
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0] == '#' || !lineStr.find("chr\t")) continue;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    T bed(v);
    vbed.push_back(bed);
  }
  return vbed;
}

template <class T>
void printBed(const vector<T> &vbed)
{
  for (auto x: vbed) {
    x.print();
    cout << endl;
  }
  cout << "bed num: " << vbed.size() << endl;
  return;
}

#endif  // READGENE_H
