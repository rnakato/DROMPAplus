/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_HAMMING_H_
#define _PW_HAMMING_H_

#include "pw_gv.h"

class shiftDist{
  double bk;
  double r;
  
 public:
  map<int, double> mp;
  map<int, double> nc;
  double nsc;
  int nsci;
  
 shiftDist(): bk(0), r(0), nsc(0), nsci(0) {}

  double getmpsum() {
    double sum(0);
    for(auto itr = mp.begin(); itr != mp.end(); ++itr) sum += itr->second;
    return sum;
  }
  
  void setControlRatio() {
    int n(0);
    for(auto itr = nc.begin(); itr != nc.end(); ++itr) {
      bk += itr->second;
      ++n;
    }
    bk /= n;
    r = 1/bk;
  }
  void getflen(int lenF3) {
    for(auto itr = mp.begin(); itr != mp.end(); ++itr) {
      if(itr->first > lenF3*2 && nsc < itr->second*r) {
	nsc = itr->second*r;
	nsci = itr->first;
      }
    }
  }
  void outputmp(ofstream &out, const Mapfile &p, const string str) {
    setControlRatio();
    getflen(p.dist.lenF3);
    double sum(getmpsum());

    out << "NSC\t"<< nsc << endl;
    out << "Estimated fragment length\t" << nsci << endl;
    out << "Background enrichment\t" << (bk*NUM_10M/p.genome.bothnread_nonred()) << endl;
    out << "Strand shift\t" << str << "\tprop\tper 10M reads\tper control" << endl;
    for(auto itr = mp.begin(); itr != mp.end(); ++itr) {
      out << itr->first << "\t" << itr->second << "\t" << (itr->second/sum)<< "\t" << (itr->second*NUM_10M/p.genome.bothnread_nonred()) << "\t" << (itr->second*r) << endl;
    }
  }
};

void hammingDist(Mapfile &p, int numthreads);
void pw_ccp(Mapfile &p, int numthreads);
void pw_Jaccard(Mapfile &p, int numthreads);

#endif /* _PW_HAMMING_H_ */
