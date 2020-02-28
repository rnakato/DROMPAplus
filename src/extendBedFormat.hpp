/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _EXTENDBEDFORMAT_HPP_
#define _EXTENDBEDFORMAT_HPP_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "../submodules/SSP/common/inline.hpp"
#include "../submodules/SSP/common/util.hpp"
#include "../submodules/SSP/common/BedFormat.hpp"

bool isStr(std::string, std::string);

class Interaction;

class GenomePosition {
 public:
  std::string chr;
  int32_t start;
 GenomePosition(): start(0) {}
  virtual ~GenomePosition(){}
  GenomePosition(const std::string &c, const std::string &s):
    chr(rmchr(c)), start(stoi(s))
  {}
  bool operator<(const GenomePosition &another) const
  {
    if (compare_chr(chr, another.chr) < 0) return 1;
    else if (compare_chr(chr, another.chr) == 0 && start < another.start) return 1;
    else return 0;
  };
};

class bed {
 public:
  std::string chr;
  int32_t start;
  int32_t end;
  int32_t summit;
 bed(): start(0), end(0), summit(0) {}
  virtual ~bed(){}
  bed(const std::string &c, const int32_t s, const int32_t e, const int32_t _summit=0):
    chr(rmchr(c)), start(s), end(e)
  {
    if (_summit) summit = _summit;
    else summit = (start + end)/2;
  }
  explicit bed(const std::vector<std::string> &s) {
    if(s.size() < 3) {
      std::cerr << "\nWarning: Bed site size < 3." << std::endl;
      return;
    }
    chr = rmchr(s[0]);
    start = stoi(s[1]);
    end = stoi(s[2]);
    summit = (start + end)/2;
  }
  void print() const { std::cout << "chr" << chr << "\t" << start  << "\t" << end ; }
  void printHead() const { std::cout << "chromosome\tstart\tend"; }
  int32_t length() const { return abs(end - start); }
  std::string getSiteStr() const {
    return "chr" + chr + "-" + std::to_string(start) + "-" + std::to_string(end);
  }

};

class bed12 : public bed {
 public:
  std::string name;
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
	  blockCount(0), blockSizes(0), blockStarts(0)
  {}
  explicit bed12(std::vector<std::string> &s):
   bed(s), rgb_r(-1), rgb_g(-1), rgb_b(-1)
  {
   int32_t num = s.size();
   if(num > 3)  name        = s[3];
   if(num > 4)  score       = stoi(s[4]);
   if(num > 5)  strand      = s[5];
   if(num > 6)  thickStart  = stoi(s[6]);
   if(num > 7)  thickEnd    = stoi(s[7]);
   if(num > 8) {
     std::vector<std::string> v;
     ParseLine(v, s[8], ',');
     if(v.size() >= 3) {
       rgb_r = stoi(v[0]);
       rgb_g = stoi(v[1]);
       rgb_b = stoi(v[2]);
     }
   }
   if(num > 9)  blockCount  = stoi(s[9]);
   if(num > 10) blockSizes  = stoi(s[10]);
   if(num > 11) blockStarts = stoi(s[11]);
 }
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
  std::string name;

 macsxls(): bed(), len(0), pileup(0), p(0), enrich(0), q(0) {}
 explicit macsxls(std::vector<std::string> &s): bed(s) {
   len    = stoi(s[3]);
   summit = stoi(s[4]);
   pileup = stod(s[5]);
   p      = stod(s[6]);
   enrich = stod(s[7]);
   q      = stod(s[8]);
   name   = s[9];
 }
 void print() const {
   std::cout << chr << "\t" << start << "\t" << end << "\t"
	<< len << "\t" << summit << "\t" << pileup << "\t"
	<< p << "\t" << enrich << "\t" << q << "\t" << name;
 }
 void printHead () const {
   std::cout << "chromosome\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname";
 }
};

class Peak : public bed {
 public:
  double pileup;
  double enrich;
  double p_inter, p_enr;
  double q;

  Peak(){}
  Peak(const std::string &c, const int32_t s, const int32_t e, const double val, const double p):
    bed(c,s,e), pileup(val), enrich(0), p_inter(p), p_enr(0), q(0) {}
  void renew(const int32_t i, const double val, const double p) {
    end = i;
    pileup += val;
    if(p_inter > p) {
      p_inter = p;
      summit = i;
    }
  }
  void print(std::ofstream &out, const int32_t id, const int32_t binsize) const {
    out << chr << "\t" << start*binsize << "\t" << end*binsize << "\t"
	<< ((end - start +1)*binsize-1) << "\t"
	<< (summit*binsize -binsize/2) << "\t" << pileup << "\t"
	<< p_inter << "\t" << enrich << "\t" << q << "\tpeak " << id << std::endl;
  }
  void printHead (std::ofstream &out) const {
    out << "chromosome\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname" << std::endl;
  }
};

template <class T=bed>
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

  void setAsMango(const std::string &lineStr) {
    if (isStr(lineStr, "color")) {
      std::cerr << "Warning: Interaction does not look mango format." << std::endl;
      return;
    }
    if (isStr(lineStr, "chrom1")) return;

    std::vector<std::string> v;
    ParseLine(v, lineStr, '\t');
    if(v.size() < 8) {
      std::cerr << "Warning: " << lineStr << " does not contain 8 columns." << std::endl;
      return;
    }
    try {
      double val, val_tmp(0);
      if(v.size() > 8) val_tmp = stod(v[15]); else val_tmp = stod(v[7]); // P
      if(val) val = -log10(val_tmp); else val = -log10(1e-12);

      vinter.emplace_back(bed({v[0], v[1], v[2]}),
			  bed({v[3], v[4], v[5]}),
			  val);
      maxval = std::max(val, maxval);
    } catch (std::exception &e) {
      PRINTERR_AND_EXIT(e.what());
    }
  }
  void setAsHICCUPS(const std::string &lineStr) {
    if (isStr(lineStr, "color")) return;
    std::vector<std::string> v;
    ParseLine(v, lineStr, '\t');
      //    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    if(v.size() < 19) {
      std::cerr << "Warning: " << lineStr << " does not contain 8 columns." << std::endl;
      return;
    }

    try {
      double val(-log10(stod(v[17]))); // fdrDonut
      vinter.emplace_back(bed({v[0], v[1], v[2]}),
			  bed({v[3], v[4], v[5]}),
			  val);
      if(std::isfinite(val)) maxval = std::max(val, maxval);
    } catch (std::exception &e) {
      PRINTERR_AND_EXIT(e.what());
    }
  }

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
	       << label << "\n"
	       << tool << "\n");

    while (!in.eof()) {
      std::string lineStr;
      getline(in, lineStr);
//      std::cout << lineStr << std::endl;
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

  void compare_bed_loop(const std::vector<bed> &bed1, const std::vector<bed> &bed2, const bool nobs) {
    int32_t aa(0), bb(0), ab(0), an(0), bn(0), nn(0), hit_bed1(0), hit_bed2(0);

    for (auto &x: bed1) {
      if (isoverlap_asBed(x, vinter)) ++hit_bed1;
    }
    for (auto &x: bed2) {
      if (isoverlap_asBed(x, vinter)) ++hit_bed2;
    }
    for (auto &x: vinter) {
      int32_t on(0);
      x.ofirst.peakovrlpd1  = isoverlap_asloop(x.first, bed1);
      x.ofirst.peakovrlpd2  = isoverlap_asloop(x.first, bed2);
      x.osecond.peakovrlpd1 = isoverlap_asloop(x.second, bed1);
      x.osecond.peakovrlpd2 = isoverlap_asloop(x.second, bed2);

      if ((x.ofirst.peakovrlpd1 && x.osecond.peakovrlpd2) || (x.ofirst.peakovrlpd2 && x.osecond.peakovrlpd1)) {
	++ab;
	++on;
      }
      if (x.ofirst.peakovrlpd1 && x.osecond.peakovrlpd1) {
	++aa;
	++on;
      }
      if (x.ofirst.peakovrlpd2 && x.osecond.peakovrlpd2) {
	++bb;
	++on;
      }
      if (on) continue;
      if (x.ofirst.peakovrlpd1 || x.osecond.peakovrlpd1) ++an;
      else if (x.ofirst.peakovrlpd2 || x.osecond.peakovrlpd2) ++bn;
      else ++nn;

    }

    std::cout << "# Number: " << vinter.size() << std::endl;
    printf("# bed1-bed2: %d (%.1f%%)\n", ab, getpercent(ab, vinter.size()));
    printf("# bed1-bed1: %d (%.1f%%)\n", aa, getpercent(aa, vinter.size()));
    printf("# bed2-bed2: %d (%.1f%%)\n", bb, getpercent(bb, vinter.size()));
    printf("# bed1-none: %d (%.1f%%)\n", an, getpercent(an, vinter.size()));
    printf("# bed2-none: %d (%.1f%%)\n", bn, getpercent(bn, vinter.size()));
    printf("# none: %d (%.1f%%)\n",      nn, getpercent(nn, vinter.size()));
    printf("# bed1: %d (%.1f%%)\n", hit_bed1, getpercent(hit_bed1, bed1.size()));
    printf("# bed2: %d (%.1f%%)\n", hit_bed2, getpercent(hit_bed2, bed2.size()));

    if(!nobs) {
      for (auto &x: vinter) {
	if((x.ofirst.peakovrlpd1 && x.osecond.peakovrlpd2) || (x.ofirst.peakovrlpd2 && x.osecond.peakovrlpd1)) x.print();
      }
    }
  }
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
    chr = rmchr(s[0]);
    start = stoi(s[1]);
    end = stoi(s[2]);
    name = s[3];
    stain = s[4];
//    std::cout << name << "," << stain << "," << start << "," << end << std::endl;
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
      //    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
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
