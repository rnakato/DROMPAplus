/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _EXTENDBEDFORMAT_HPP_
#define _EXTENDBEDFORMAT_HPP_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "../../submodules/SSP/common/inline.hpp"
#include "../../submodules/SSP/common/util.hpp"

class bed {
public:
  std::string chr;
  int32_t start;
  int32_t end;
  int32_t summit;
  std::string name;

  bed(): start(0), end(0), summit(0) {}
  virtual ~bed(){}

  bed(const std::string &c, const int32_t s, const int32_t e, const int32_t _summit=0):
    chr(rmchr(c)), start(s), end(e), name("")
  {
    if (_summit) summit = _summit;
    else summit = (start + end)/2;
  }

  explicit bed(const std::vector<std::string> &s) {
    if(s.size() < 3) {
      std::cerr << "\nWarning: Bed site size < 3." << std::endl;
      return;
    }

    try {
      chr = rmchr(s[0]);
      start = stoi(s[1]);
      end = stoi(s[2]);
      summit = (start + end)/2;
      if(s.size() >= 4) name = s[3];
    } catch (std::exception &e) {
      PRINTERR_AND_EXIT("invalid columns in BED format. " + std::string(e.what()));
    }

  }
  void print() const { std::cout << "chr" << chr << "\t" << start  << "\t" << end ; }
  void printHead() const { std::cout << "chromosome\tstart\tend"; }
  int32_t length() const { return abs(end - start); }
  std::string getSiteStr() const {
    return "chr" + chr + "-" + std::to_string(start) + "-" + std::to_string(end);
  }
  std::string getSiteStrTAB() const {
    return "chr" + chr + "\t" + std::to_string(start) + "\t" + std::to_string(end);
  }
  std::string getSiteStrTABwithNAME() const {
    return getSiteStrTAB() + "\t" + name;
  }
};


class GenomicPosition {
public:
  std::string chr;
  int32_t start;

  GenomicPosition(): start(0) {}
  GenomicPosition(const std::string &c, const std::string &s):
    chr(rmchr(c)), start(stoi(s))
  {}

  bool operator<(const GenomicPosition &another) const
  {
    if (compare_chr(chr, another.chr) < 0) return 1;
    else if (compare_chr(chr, another.chr) == 0 && start < another.start) return 1;
    else return 0;
  };
};

class bed12 : public bed {
public:
//  std::string name;
  int32_t score;
  std::string strand;
  int32_t thickStart;
  int32_t thickEnd;
  int32_t rgb_r, rgb_g, rgb_b;
  int32_t blockCount;
  int32_t blockSizes;
  int32_t blockStarts;
  bed12(): bed(), score(0), thickStart(0), thickEnd(0),
           rgb_r(0), rgb_g(0), rgb_b(0),
           blockCount(0), blockSizes(0), blockStarts(0) {}
  bed12(std::vector<std::string> &s);

  void print() const {
    std::cout << chr << "\t" << start << "\t" << end << "\t"
              << name << "\t" << score << "\t" << strand << "\t"
              << thickStart << "\t" << thickEnd << "\t"
              << rgb_r << "\t" << rgb_g << "\t" << rgb_b << "\t"
              << blockCount << "\t" << blockSizes << "\t" << blockStarts;
  }
  void printHead() const {
    std::cout << "chromosome\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\tRgb_r\tRgb_g\tRgb_b\tblockCount\tblockSizes\tblockStarts";
  }
};

class macsxls : public bed {
public:
  int32_t len;
  double pileup;
  double p;
  double enrich;
  double q;
//  std::string name;

  macsxls(): bed(), len(0), pileup(0), p(0), enrich(0), q(0) {}
  explicit macsxls(std::vector<std::string> &s): bed(s) {

    try {
      len    = stoi(s[3]);
      summit = stoi(s[4]);
      pileup = stod(s[5]);
      p      = stod(s[6]);
      enrich = stod(s[7]);
      q      = stod(s[8]);
      name   = s[9];
    } catch (std::exception &e) {
      PRINTERR_AND_EXIT("invalid columns in macsxls format. " + std::string(e.what()));
    }
  }
  void print() const {
    std::cout << chr << "\t" << start << "\t" << end << "\t"
              << len << "\t" << summit << "\t" << pileup << "\t"
              << p << "\t" << enrich << "\t" << q << "\t" << name;
  }
  void printHead () const {
    std::cout << "chromosome\tstart\tend\tlength\tabs_summit\tpileup\t log10(pvalue)\tfold_enrichment\t log10(qvalue)\tname";
  }
};

class Peak : public bed {
  int32_t binsize;

public:
  double pileup, pileup_input;
  double p_inter, p_enr;

  Peak(){}
  Peak(const std::string &c, const int32_t _binsize,
       const int32_t s, const int32_t e,
       const double val, const double _p_inter,
       const double val_input=0, const double _p_enr=0):
    bed(c, s, e), binsize(_binsize), pileup(val), pileup_input(val_input),
    p_inter(_p_inter), p_enr(_p_enr)
  {}

  void renew(const int32_t e, const double val, const double _p_inter, const double val_input=0, const double _p_enr=0) {
    end = e;
    pileup += val;
    pileup_input += val_input;

    if(p_inter < _p_inter) {
      p_inter = _p_inter;
      summit = e - binsize/2;
    }
    if(p_enr < _p_enr) p_enr = _p_enr;
  }

  double getenrichment() const { return getratio(pileup, pileup_input); }

  void print(std::ofstream &out, const int32_t id) const {
    out << "chr" << chr << "\t" << start << "\t" << end << "\t"
        << length() << "\t" << summit << "\t"
        << std::fixed << std::setprecision(2)
        << pileup << "\t" << pileup_input << "\t" << getenrichment() << "\t"
        << p_inter << "\t" << p_enr << "\tpeak " << id << std::endl;
  }

  void print() const {
    std::cout << "chr"<< chr << "\t" << start << "\t" << end << "\t"
              << length() << "\t"
              << summit << "\t"
              << std::fixed << std::setprecision(2)
              << pileup << "\t" << pileup_input << "\t" << getenrichment() << "\t"
              << p_inter << "\t" << p_enr << "\tpeak " << std::endl;
  }

  void printHead (std::ofstream &out) const {
    out << "chromosome\tstart\tend\tlength\tsummit\t"
        << "pileup (ChIP)\tpileup (Input)\tEnrichment\t"
        << "-log10(p_internal)\t-log10(p_enrichment)\tname" << std::endl;
  }
};

template <class T>
class vbed {
  std::vector<T> bed;
  std::string label;

public:
  vbed(){}
  vbed(const std::vector<T> &v, const std::string &l):
    bed(v), label(l)
  {}
  const std::string & getlabel() const { return label; }
  std::vector<T> getvBed() const { return bed; }
};


template <class T>
std::vector<T> parseBed(const std::string &fileName)
{
  DEBUGprint_FUNCStart();

  std::vector<T> vbed;
  std::ifstream in(fileName);
  if(!in) {
    PRINTERR_AND_EXIT("Error: BED file " << fileName << " does not exist.");
  }

  while (!in.eof()) {
    std::string lineStr;
    getline(in, lineStr);

    if(lineStr.empty() || lineStr[0] == '#') continue;
    if(lineStr.find("chromosome") != std::string::npos) continue;  // header
    if(lineStr.find("description") != std::string::npos) continue;  // ChromHMM.dense.bed

    std::vector<std::string> v;
    ParseLine(v, lineStr, '\t');

    if(v[1] == "start") continue;
    vbed.emplace_back(v);
  }

  DEBUGprint_FUNCend();
  return vbed;
}

template <class T>
void printBed(const std::vector<T> &vbed)
{
  for (auto &x: vbed) {
    x.print();
    std::cout << std::endl;
  }
  std::cout << "Total number: " << vbed.size() << std::endl;
  return;
}

class peakoverlapped {
public:
  bool peakovrlpd1;
  bool peakovrlpd2;
  peakoverlapped(): peakovrlpd1(false), peakovrlpd2(false){}
};

class Interaction {
  double val;
public:
  bed first;
  bed second;
  peakoverlapped ofirst;
  peakoverlapped osecond;
  Interaction(): val(0) {}
  Interaction(const std::string &c1, const int32_t s1, const int32_t e1,
              const std::string &c2, const int32_t s2, const int32_t e2,
              const double v=0):
    val(v), first(c1, s1, e1), second(c2, s2, e2)
  {}
  Interaction(const bed &b1, const bed &b2, const double _val):
    val(_val), first(b1), second(b2)
  {}

  double getval() const { return val; }
  void print() const {
    first.print();
    std::cout << "\t";
    second.print();
    std::cout << "\t" << val << std::endl;
  }
};

class InteractionSet {
  std::vector<Interaction> vinter;
  double maxval;
  std::string label;

  void setAsMango(const std::string &lineStr);
  void setAsHICCUPS(const std::string &lineStr);

public:
  InteractionSet(const std::string &fileName, const std::string &l, const std::string &tool):
    maxval(1e-10), label(l)
  {
    std::ifstream in(fileName);
    if(!in) {
      PRINTERR_AND_EXIT("Error: Interaction file " << fileName << " does not exist.");
    }
    DEBUGprint("Add InteractionSet.. (--inter)\n"
               << fileName << "\n"
               << label    << "\n"
               << tool     << "\n");

    while (!in.eof()) {
      std::string lineStr;
      getline(in, lineStr);
      if(lineStr.empty() || lineStr[0] == '#') continue;
      if (tool == "mango") setAsMango(lineStr);
      else setAsHICCUPS(lineStr);
    }
    //    print();
  }
  const std::vector<Interaction> & getvinter() const { return vinter; }
  const std::string & getlabel() const { return label; }
  double getmaxval() const { return maxval; }
  size_t getnum() const { return vinter.size(); }
  void print() const {
    printBed(vinter);
    std::cout << "maxval: " << maxval << std::endl;
  }

  bool isoverlap_asloop(const bed &loop, const std::vector<bed> &bed) const {
    for (auto &b: bed) {
      if (loop.chr == b.chr && my_overlap(loop.start, loop.end, b.start, b.end)) return true;
    }
    return false;
  }

  bool isoverlap_asBed(const bed &bed, const std::vector<Interaction> &vinter) const {
    for (auto &x: vinter) {
      if ((x.first.chr  == bed.chr && my_overlap(x.first.start,  x.first.end,  bed.start, bed.end)) ||
          (x.second.chr == bed.chr && my_overlap(x.second.start, x.second.end, bed.start, bed.end)))
        return true;
    }
    return false;
  }

  void compare_bed_loop(const std::vector<bed> &bed1, const std::vector<bed> &bed2, const bool nobs);

};

class cytoband {
public:
  std::string chr;
  int32_t start;
  int32_t end;
  std::string name;
  std::string stain;
  cytoband(): start(0), end(0) {}
  virtual ~cytoband(){}
  explicit cytoband(const std::vector<std::string> &s) {
    if(s.size() < 5) {
      std::cerr << "\nWarning: Cytoband size < 5." << std::endl;
      return;
    }

    try {
      chr = rmchr(s[0]);
      start = stoi(s[1]);
      end = stoi(s[2]);
      name = s[3];
      stain = s[4];
      //    std::cout << name << "," << stain << "," << start << "," << end << std::endl;
    } catch (std::exception &e) {
      PRINTERR_AND_EXIT("invalid columns in cytoband format. " + std::string(e.what()));
    }

  }

  void print() const {
    std::cout << "chr" << chr << "\t" << start  << "\t" << end
              << "\t" << name << "\t" << stain << std::endl;
  }
  int32_t getlen() const { return end - start; }
};


template <class T>
std::unordered_map<std::string, std::vector<T>> parseBed_Hash(const std::string &fileName)
{
  std::unordered_map<std::string, std::vector<T>> bedmap;
  std::ifstream in(fileName);
  if(!in) PRINTERR_AND_EXIT("Error: BED file does not exist.");

  std::string lineStr;
  std::vector<std::string> v;
  while (!in.eof()) {
    getline(in, lineStr);

    if(lineStr.empty() || lineStr[0] == '#') continue;
    ParseLine(v, lineStr, '\t');
    if(v[1] == "start") continue;
    bedmap[rmchr(v[0])].emplace_back(bed(v));
  }

  return bedmap;
}

template <class T>
void printBed_Hash(const std::unordered_map<std::string, std::vector<T>> &mp)
{
  for(auto vbed: mp) printBed(vbed.second);
  return;
}

#endif  // _EXTENDBEDFORMAT_HPP_
