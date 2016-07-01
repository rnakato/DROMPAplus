/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_SHIFTPROFILE_P_H_
#define _PW_SHIFTPROFILE_P_H_

#include "pw_gv.h"

#define HD_FROM 200
#define HD_TO   600
#define NG_FROM 4000
#define NG_TO   5000
#define NG_STEP 100

vector<char> genVector(const strandData &seq, int start, int end);
boost::dynamic_bitset<> genBitset(const strandData &seq, int start, int end, long chrlen);

class _shiftDist {
  int lenF3;
 public:
  map<int, double> mp;
  map<int, double> nc;
  long nread;
  double r;
  double bk;
  double nsc;
  int nsci;
  double rchr;
  
 _shiftDist(const Mapfile &p, long n=0): lenF3(p.dist.lenF3), nread(n), r(0), bk(0), nsc(0), nsci(0), rchr(1) {}
  
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
  void setflen() {
    for(auto itr = mp.begin(); itr != mp.end(); ++itr) {
      if(itr->first > lenF3*2 && nsc < itr->second*r) {
	nsc = itr->second*r;
	nsci = itr->first;
      }
    }
  }
  double getBackEnrich(long nread) {
    return bk *NUM_10M /(double)nread;
  }
  void outputmp(const string filename, long nread, string name) {
    setControlRatio();
    setflen();
    double sum(getmpsum());
    
    ofstream out(filename);
    out << "NSC\t"<< nsc << endl;
    out << "Estimated fragment length\t" << nsci << endl;
    out << "Background enrichment\t" << getBackEnrich(nread) << endl;
    out << "Strand shift\t" << name << "\tprop\tper 10M reads\tper control" << endl;
    for(auto itr = mp.begin(); itr != mp.end(); ++itr) {
      out << itr->first << "\t" << itr->second << "\t" << (itr->second/sum)<< "\t" << (itr->second*NUM_10M/(double)nread) << "\t" << (itr->second * r) << endl;
    }
  }
};

class shiftDist {
 protected:
  int mp_from;
  int mp_to;
  int mp_step;
  int ng_from;
  int ng_to;
  int ng_step;
  
 public:
  string name;
  _shiftDist genome;
  vector<_shiftDist> chr;
  
 shiftDist(string n, const Mapfile &p): mp_from(HD_FROM), mp_to(HD_TO), mp_step(1), ng_from(NG_FROM), ng_to(NG_TO), ng_step(NG_STEP), name(n), genome(p) {
    for(auto x:p.chr) {
      if(x.isautosome()) genome.nread += x.bothnread_nonred();
    }
    for(auto x:p.chr) {
      _shiftDist v(p, x.bothnread_nonred());
      v.rchr = v.nread/(double)genome.nread;
      cout << v.rchr << "\t"<< v.nread << "\t"<< genome.nread << endl;
      chr.push_back(v);
    }
  }
  void add2genome(const _shiftDist &x, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
#ifdef DEBUG
    cout << "add2genome" << endl;
#endif
    addmp(genome.mp, x.mp, x.rchr);
    addmp(genome.nc, x.nc, x.rchr);
  }
};

class shiftJacVec : public shiftDist {
 public:
 shiftJacVec(const Mapfile &p): shiftDist("Jaccard index", p) {}
  
  void setDist(int i, const vector<char> &fwd, const vector<char> &rev, int start, int end);
  void execchr(const Mapfile &p, int i, int start, int end, boost::mutex &mtx) {
    auto fwd = genVector(p.chr[i].seq[STRAND_PLUS],  start, end);
    auto rev = genVector(p.chr[i].seq[STRAND_MINUS], start, end);
    
    setDist(i, fwd, rev, start, end);
  }
};

class shiftJacBit : public shiftDist {
 public:
 shiftJacBit(const Mapfile &p): shiftDist("Jaccard index", p) {}
  
  void setDist(int i, boost::dynamic_bitset<> &fwd, const boost::dynamic_bitset<> &rev, int start, int end);
  void execchr(const Mapfile &p, int i, int start, int end, boost::mutex &mtx) {
    auto fwd = genBitset(p.chr[i].seq[STRAND_PLUS],  start, end, p.chr[i].len);
    auto rev = genBitset(p.chr[i].seq[STRAND_MINUS], start, end, p.chr[i].len);
    
    setDist(i, fwd, rev, start, end);
  }
};

class shiftCcp : public shiftDist {
 public:
 shiftCcp(const Mapfile &p): shiftDist("Cross correlation", p) {}
  
  void setDist(int i, const vector<char> &fwd, const vector<char> &rev, int start, int end);
  void execchr(const Mapfile &p, int i, int start, int end, boost::mutex &mtx) {
    auto fwd = genVector(p.chr[i].seq[STRAND_PLUS],  start, end);
    auto rev = genVector(p.chr[i].seq[STRAND_MINUS], start, end);
    
    setDist(i, fwd, rev, start, end);
  }
};

class shiftHamming : public shiftDist {
 public:
 shiftHamming(const Mapfile &p): shiftDist("Hamming distance", p) {}

  void setDist(int i, boost::dynamic_bitset<> &fwd, const boost::dynamic_bitset<> &rev, int start, int end);
  void execchr(const Mapfile &p, int i, int start, int end, boost::mutex &mtx) {
    auto fwd = genBitset(p.chr[i].seq[STRAND_PLUS],  start, end, p.chr[i].len);
    auto rev = genBitset(p.chr[i].seq[STRAND_MINUS], start, end, p.chr[i].len);
    
    setDist(i, fwd, rev, start, end);
  }
};

#endif /* _PW_SHIFTPROFILE_P_H_ */
