/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _COMMON_H_
#define _COMMON_H_

#include <iostream>
#include <string>
#include "seq.h"

#define VERSION "3.3.0+"
#define NUM_1K 1000
#define NUM_1M 1000000
#define NUM_10M 10000000
#define NUM_100M 100000000

#define VALUE2WIGARRAY(v) ((v) * 1000.0)
#define WIGARRAY2VALUE(v) ((v) * (1.0/1000.0))

vector<string> str_wigfiletype = {"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};


enum PWfile_Type {
  TYPE_BINARY,
  TYPE_COMPRESSWIG,
  TYPE_UNCOMPRESSWIG,
  TYPE_BEDGRAPH,
  TYPE_BIGWIG,
  PWFILETYPENUM
};

class Peak : public bed {
 public:
  int summit;
  double pileup;
  double enrich;
  double p_inter, p_enr;
  double q;
 Peak(int s, int e, string c, double val, double p):
  bed(s,e,c), summit(s), pileup(val), enrich(0), p_inter(p), p_enr(0), q(0) {}
  void renew(int i, double val, double p) {
    end = i;
    pileup += val;
    if(p_inter > p) {
      p_inter = p;
      summit = i;
    }
  }
  void print(ofstream &out, int id, int binsize) const {
    out << chr << "\t" << start*binsize << "\t" << end*binsize << "\t"
	<< ((end - start +1)*binsize-1) << "\t" << (summit*binsize -binsize/2) << "\t" << pileup << "\t"
	<< p_inter << "\t" << enrich << "\t" << q << "\tpeak " << id << endl;
  }
  void printHead (ofstream &out) const {
    out << "chromosome\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname" << endl;
  }
};

template <class T>
class var {
  std::string name;
 public:
  T val;
 var(std::string str, T v=0): name(str), val(v) {}
  friend std::ostream &operator<<(std::ostream &stream, var v) {
    return stream << v.name << "\t" << v.val;
  }
  var operator++() {
    val++;
    return *this;
  }
};

template <class T>
class temp {
  int on;
  std::string str;
 temp(): on(0){}
  void print() {
    if(!on) {
      std::cout << str << std::endl;
      on++;
    }
  }
};

#endif /* _COMMON_H_ */
