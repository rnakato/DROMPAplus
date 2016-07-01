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
boost::dynamic_bitset<> genBitset(const strandData &seq, int, int);

class _shiftDist {
  int lenF3;
 public:
  map<int, double> mp;
  map<int, double> nc;
  int start;
  int end;
  int width4mp;
  long nread;
  double r;
  double bk;
  double nsc;
  int nsci;
  double rchr;

  int numthreads;
  
 _shiftDist(const Mapfile &p, int s=0, int e=0, int ng_to=0, long n=0, int nthre=0): lenF3(p.dist.lenF3), start(s), end(e), width4mp(e-s-ng_to), nread(n), r(0), bk(0), nsc(0), nsci(0), rchr(1), numthreads(nthre) {}

  void setmp(int i, double val, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    mp[i] = val;
  }
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
  int ng_from;
  int ng_to;
  int ng_step;
  vector<range> seprange;
  
 public:
  string name;
  _shiftDist genome;
  vector<_shiftDist> chr;
  
 shiftDist(string n, const Mapfile &p, int numthreads): mp_from(HD_FROM), mp_to(HD_TO), ng_from(NG_FROM), ng_to(NG_TO), ng_step(NG_STEP), name(n), genome(p) {
    for(auto x:p.chr) {
      if(x.isautosome()) genome.nread += x.bothnread_nonred();
    }
    for(auto x:p.chr) {
      _shiftDist v(p, 0, x.len, ng_to, x.bothnread_nonred(), numthreads);
      v.rchr = v.nread/(double)genome.nread;
      chr.push_back(v);
    }
    // seprange
    int length(mp_to+mp_from);
    int sepsize = length/numthreads +1;
    for(int i=0; i<numthreads; i++) {
      int s = i*sepsize;
      int e = (i+1)*sepsize;
      if(i==numthreads-1) e = length;
      range sep(s - mp_from, e - mp_from);
      seprange.push_back(sep);
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
 shiftJacVec(const Mapfile &p, int numthreads): shiftDist("Jaccard index", p, numthreads) {}
  
  void setDist(_shiftDist &chr, const vector<char> &fwd, const vector<char> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genVector(p.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genVector(p.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftJacBit : public shiftDist {
 public:
 shiftJacBit(const Mapfile &p, int numthreads): shiftDist("Jaccard index", p, numthreads) {}

  void setDist(_shiftDist &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genBitset(p.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genBitset(p.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftCcp : public shiftDist {
 public:
 shiftCcp(const Mapfile &p, int numthreads): shiftDist("Cross correlation", p, numthreads) {}
  
  void setDist(_shiftDist &chr, const vector<char> &fwd, const vector<char> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genVector(p.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genVector(p.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftHamming : public shiftDist {
 public:
 shiftHamming(const Mapfile &p, int numthreads): shiftDist("Hamming distance", p, numthreads) {}

  void setDist(_shiftDist &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genBitset(p.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genBitset(p.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);
    
    setDist(chr[i], fwd, rev);
  }
};

#endif /* _PW_SHIFTPROFILE_P_H_ */
