/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _BEDFORMAT_HPP_
#define _BEDFORMAT_HPP_

#include <fstream>
#include <boost/algorithm/string.hpp>

std::string rmchr(const std::string &chr);

class bed {
 public:
  std::string chr;
  int32_t start;
  int32_t end;
  int32_t summit;
 bed(): start(0), end(0), summit(0) {}
  virtual ~bed(){}
  bed(const std::string &c, const int32_t s, const int32_t e):
    chr(rmchr(c)), start(s), end(e) {}
  bed(const std::vector<std::string> &s):
    chr(rmchr(s[0])), start(stoi(s[1])), end(stoi(s[2])), summit((start + end)/2) {}
  void print()     const { std::cout << "chr" << chr << "\t" << start  << "\t" << end ; }
  void printHead() const { std::cout << "chromosome\tstart\tend"; }
  int32_t length() const { return abs(end - start); }
};

class bed12 : public bed {
 public:
  std::string name;
  int32_t score;
  std::string strand;
  int32_t thickStart;
  int32_t thickEnd;
  std::string itemRgb;
  int32_t blockCount;
  int32_t blockSizes;
  int32_t blockStarts;
 bed12(): bed() {}
 bed12(std::vector<std::string> s): bed(s) {
   int32_t num = s.size();
   if(num > 3)  name        = s[3];
   if(num > 4)  score       = stoi(s[4]);
   if(num > 5)  strand      = s[5];
   if(num > 6)  thickStart  = stoi(s[6]);
   if(num > 7)  thickEnd    = stoi(s[7]);
   if(num > 8)  itemRgb     = s[8];
   if(num > 9)  blockCount  = stoi(s[9]);
   if(num > 10) blockSizes  = stoi(s[10]);
   if(num > 11) blockStarts = stoi(s[11]);
 }
 void print() const {
   std::cout << chr << "\t" << start << "\t" << end << "\t"
	<< name << "\t" << score << "\t" << strand << "\t"
	<< thickStart << "\t" << thickEnd << "\t" << itemRgb << "\t"
	<< blockCount << "\t" << blockSizes << "\t" << blockStarts;
 }
 void printHead() const {
   std::cout << "chromosome\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts";
 }
};

class macsxls : public bed {
 public:
  int32_t len;
  int32_t summit;
  double pileup;
  double p;
  double enrich;
  double q;
  std::string name;

 macsxls(): bed() {}
 macsxls(std::vector<std::string> s): bed(s) {
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
  int32_t summit;
  double pileup;
  double enrich;
  double p_inter, p_enr;
  double q;
 Peak(const std::string &c, const int32_t s, const int32_t e, const double val, const double p):
   bed(c,s,e), summit(s), pileup(val), enrich(0), p_inter(p), p_enr(0), q(0) {}
  void renew(int32_t i, double val, double p) {
    end = i;
    pileup += val;
    if(p_inter > p) {
      p_inter = p;
      summit = i;
    }
  }
  void print(std::ofstream &out, int32_t id, int32_t binsize) const {
    out << chr << "\t" << start*binsize << "\t" << end*binsize << "\t"
	<< ((end - start +1)*binsize-1) << "\t"
	<< (summit*binsize -binsize/2) << "\t" << pileup << "\t"
	<< p_inter << "\t" << enrich << "\t" << q << "\tpeak " << id << std::endl;
  }
  void printHead (std::ofstream &out) const {
    out << "chromosome\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname" << std::endl;
  }
};

template <class T>
std::vector<T> parseBed(const std::string &fileName)
{
  std::vector<T> vbed;
  std::ifstream in(fileName);
  if(!in) {
    std::cerr << "Error: BED file does not exist." << std::endl;
    std::exit(1);
  }

  std::string lineStr;
  std::vector<std::string> v;
  while (!in.eof()) {
    getline(in, lineStr);
    
    if(lineStr.empty() || lineStr[0] == '#') continue;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    if(v[1] == "start") continue;
    vbed.emplace_back(v);
  }

  return vbed;
}

template <class T>
void printBed(const std::vector<T> &vbed)
{
  for (auto &x: vbed) {
    x.print();
    std::cout << std::endl;
  }
  std::cout << "bed num: " << vbed.size() << std::endl;
  return;
}

#endif  // _BEDFORMAT_HPP_
