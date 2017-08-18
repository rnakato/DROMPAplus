/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _SEQSTATS_HPP_
#define _SEQSTATS_HPP_

#include <fstream>
#include <boost/thread.hpp>
#include <boost/algorithm/string.hpp>
#include "../common/seq.hpp"
#include "../common/inline.hpp"
#include "ReadBpStatus.hpp"

template <class Thead, class... Tbody>
void printList(Thead head, Tbody... body);

  
class strandData {
 public:
  std::vector<Read> vRead;
  uint64_t nread_nonred;
  uint64_t nread_red;
  double nread_rpm;
  double nread_afterGC;

 strandData(): nread_nonred(0), nread_red(0), nread_rpm(0), nread_afterGC(0) {}
  size_t getnread() const { return vRead.size(); }
  void print() const {
    printList(getnread(), nread_nonred, nread_red, nread_rpm, nread_afterGC);
  }

  void setnread_nonread_nofilter() {
    nread_nonred = getnread();
  }
};

template <class T>
void printSeqStats(const T &obj)
{
  printList(obj.getname(), obj.getlen(), obj.getlenmpbl(),
	    obj.getnread(Strand::BOTH),
	    obj.getnread_nonred(Strand::BOTH),
	    obj.getnread_red(Strand::BOTH),
	    obj.getnread_rpm(Strand::BOTH),
	    obj.getnread_afterGC(Strand::BOTH),
	    obj.getdepth());
}

class SeqStats {
  enum {STRANDNUM=2};

  std::string name;
  uint64_t len, len_mpbl;
  strandData seq[STRANDNUM];
  bool Greekchr;
  double depth;
  uint64_t nread_inbed;
  double sizefactor;

 public:
  
 SeqStats(std::string s, int32_t l):
  name(rmchr(s)), len(l), len_mpbl(l), 
  Greekchr(false), depth(0), nread_inbed(0),
  sizefactor(0) {}

  const std::string & getname() const { return name; }
  void addfrag(const Fragment &frag) {
    seq[frag.strand].vRead.emplace_back(frag, len);
  }
  uint64_t getnread (const Strand::Strand strand) const {
    if(strand==Strand::BOTH) return seq[Strand::FWD].getnread() + seq[Strand::REV].getnread();
    else return seq[strand].getnread();
  }
  uint64_t getnread_nonred (const Strand::Strand strand) const {
    if(strand==Strand::BOTH) return seq[Strand::FWD].nread_nonred + seq[Strand::REV].nread_nonred;
    else return seq[strand].nread_nonred;
  }
  uint64_t getnread_red (const Strand::Strand strand) const {
    if(strand==Strand::BOTH) return seq[Strand::FWD].nread_red + seq[Strand::REV].nread_red;
    else return seq[strand].nread_red;
  }
  uint64_t getnread_rpm (const Strand::Strand strand) const {
    if(strand==Strand::BOTH) return seq[Strand::FWD].nread_rpm + seq[Strand::REV].nread_rpm;
    else return seq[strand].nread_rpm;
  }
  uint64_t getnread_afterGC (const Strand::Strand strand) const {
    if(strand==Strand::BOTH) return seq[Strand::FWD].nread_afterGC + seq[Strand::REV].nread_afterGC;
    else return seq[strand].nread_afterGC;
  }
  uint64_t getnread_inbed() const { return nread_inbed; }

  void setnread_nonread_nofilter() {
    for (auto strand: {Strand::FWD, Strand::REV}) seq[strand].setnread_nonread_nofilter();
  }
  
  void incNReadNonred (const Strand::Strand strand) { ++seq[strand].nread_nonred; }
  void incNReadRed    (const Strand::Strand strand) { ++seq[strand].nread_red;    }

  const std::vector<Read> & getvReadref (const Strand::Strand strand) const {
    return seq[strand].vRead;
  }
  std::vector<Read> & getvReadref_notconst (const Strand::Strand strand) {
    return seq[strand].vRead;
  }
  void setdepth(const double d) { depth = d; }
  void addReadAfterGC(const Strand::Strand strand, const double w, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    seq[strand].nread_afterGC += w;
  }
  
  uint64_t getlen()       const { return len; }
  uint64_t getlenmpbl()   const { return len_mpbl; }
  double   getpmpbl()     const { return getratio(getlenmpbl(), getlen()); }
  double   getsizefactor()const { return sizefactor; }
  double   getdepth()     const { return depth; }

  void setF5ToRead(const int32_t flen) {
    int32_t d;
    for (auto strand: {Strand::FWD, Strand::REV}) {
      if(strand == Strand::FWD) d = flen; else d = -flen;
      for(auto &x: seq[strand].vRead) x.F5 = x.F3 + d;
    }
  }
  
  void printvRead() const {
#ifdef PRINTREAD
    for (auto strand: {Strand::FWD, Strand::REV}) {
      for(auto &x: seq[strand].vRead) {
	std::cout << "chr:"  << getname() << "\t"
		  << "strand:"  << strand << "\t";
	x.print();
      }
    }
#endif 
  }

  void setFRiP(const std::vector<bed> &vbed) {
    std::vector<BpStatus> array(getlen(), BpStatus::MAPPABLE);
    OverrideBedToArray(array, getname(), vbed);
    
    for (auto strand: {Strand::FWD, Strand::REV}) {
      for (auto &x: seq[strand].vRead) {
	if(x.duplicate) continue;
	int32_t s(std::min(x.F3, x.F5));
	int32_t e(std::max(x.F3, x.F5));
	//	std::cout << x.F3 << "\t" << x.F5 <<"\t" << s << "\t" << e << "\n";
	for(int32_t i=s; i<=e; ++i) {
	  if(array[i] == BpStatus::INBED) {
	    x.inpeak = 1;
	    ++nread_inbed;
	    break;
	  }
	}
      }
    }
    return;
  }

  double getFRiP() const {
    return getratio(nread_inbed, getnread_nonred(Strand::BOTH));
  }
  void setsizefactor(const double w) {
    sizefactor = w;
    for (auto strand: {Strand::FWD, Strand::REV}) seq[strand].nread_rpm = seq[strand].nread_nonred * sizefactor;
  }
  void Greekchron() { Greekchr = true; }

  bool isautosome() const {
    int32_t chrnum(0);
    try {
      chrnum = stoi(name);
    } catch (std::invalid_argument e) {  // 数値以外
      if(Greekchr) { 
	if(name=="I")         chrnum = 1;
	else if(name=="II")   chrnum = 2;
	else if(name=="III")  chrnum = 3;
	else if(name=="IV")   chrnum = 4;
	else if(name=="V")    chrnum = 5;
	else if(name=="VI")   chrnum = 6;
	else if(name=="VII")  chrnum = 7;
	else if(name=="VIII") chrnum = 8;
	else if(name=="IX")   chrnum = 9;
	else if(name=="X")    chrnum = 10;
	else if(name=="XI")   chrnum = 11;
	else if(name=="XII")  chrnum = 12;
	else if(name=="XIII") chrnum = 13;
	else if(name=="XIV")  chrnum = 14;
	else if(name=="XV")   chrnum = 15;
	else if(name=="XVI")  chrnum = 16;
      }
      if(name=="2L") chrnum = 1;
      if(name=="2R") chrnum = 2;
      if(name=="3L") chrnum = 3;
      if(name=="3R") chrnum = 4;
    }
    if(chrnum) return true;
    else       return false;
  }

  void getMptable(const std::string &mptable)
  {
    std::string lineStr;
    std::vector<std::string> v;
    std::ifstream in(mptable);
    if(!in) PRINTERR("Could nome open " << mptable << ".");
    while (!in.eof()) {
      getline(in, lineStr);
      if(lineStr.empty() || lineStr[0] == '#') continue;
      boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
      if(name == rmchr(v[0])) len_mpbl = stoi(v[1]);
    }
    return;
  }
};

#endif /* _SEQSTATS_HPP_ */
