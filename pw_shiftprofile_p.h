
/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_SHIFTPROFILE_P_H_
#define _PW_SHIFTPROFILE_P_H_

#include "pw_gv.h"

#define MP_FROM 200
#define MP_TO   600
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
  int end4mp;
  long nread;
  long len;
  double r;
  double bk;
  double nsc;
  int nsci;
  double rchr;

 _shiftDist(const Mapfile &p, int s=0, int e=0, long n=0, long l=0, double w=1): lenF3(p.dist.lenF3), start(s), end(e), end4mp(e-s-NG_TO), nread(n), len(l), r(0), bk(0), nsc(0), nsci(0), rchr(1) {}

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
  void setflen(double w) {
    int shift_min(lenF3*1.2);
    nsc = mp[shift_min+1]*w;
    for(auto itr = mp.begin(); itr != mp.end(); ++itr) {
      if(itr->first > shift_min && nsc < itr->second*r*w) {
	nsc = itr->second*r*w;
	nsci = itr->first;
      }
    }
  }
  void outputmp(const string filename, string name, double weight) {
    setControlRatio();
    setflen(weight);
    double sum(getmpsum());
    double rRPKM = (NUM_10M/(double)nread) / (NUM_100M/(double)len);
    double be(bk * rRPKM);
    double const_bu(0.024);
    
    ofstream out(filename);
    out << "NSC\t"<< nsc*weight << endl;
    out << "Estimated fragment length\t" << nsci << endl;
    out << "Background enrichment\t" << be << endl;
    out << "Background uniformity\t" << const_bu / be << endl;

    out << "Strand shift\t" << name << "\tprop\tper 10M reads\tper control" << endl;
    for(auto itr = mp.begin(); itr != mp.end(); ++itr) 
      out << itr->first            << "\t"
	  << itr->second           << "\t"
	  << (itr->second/sum)     << "\t"
	  << (itr->second * rRPKM) << "\t"
	  << (itr->second * r)     << endl;
    
  }
};

class shiftDist {
 protected:
  vector<range> seprange;
  
 public:
  string name;
  _shiftDist genome;
  vector<_shiftDist> chr;
  
 shiftDist(string n, const Mapfile &p, int numthreads): name(n), genome(p) {
    for(auto x:p.chr) {
      if(x.isautosome()) {
	genome.nread += x.bothnread_nonred();
	genome.len   += x.len_mpbl;
      }
    }
    for(auto x:p.chr) {
      _shiftDist v(p, 0, x.len, x.bothnread_nonred(), x.len_mpbl);
      v.rchr = v.nread/(double)genome.nread;
      chr.push_back(v);
    }
    // seprange
    int length(MP_TO+MP_FROM);
    int sepsize = length/numthreads +1;
    for(int i=0; i<numthreads; i++) {
      int s = i*sepsize;
      int e = (i+1)*sepsize;
      if(i==numthreads-1) e = length;
      range sep(s - MP_FROM, e - MP_FROM);
      seprange.push_back(sep);
    }
  }
  void add2genome(const _shiftDist &x, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
#ifdef DEBUG
    cout << "add2genome.." << flush;
#endif
    addmp(genome.mp, x.mp, x.rchr);
    addmp(genome.nc, x.nc, x.rchr);
  }
};

class shiftJacVec : public shiftDist {
 public:
  double w;
 shiftJacVec(const Mapfile &p, int numthreads): shiftDist("Jaccard index", p, numthreads), w(1) {}

  void setDist(_shiftDist &chr, const vector<char> &fwd, const vector<char> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genVector(p.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genVector(p.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftJacBit : public shiftDist {
 public:
  double w;
 shiftJacBit(const Mapfile &p, int numthreads): shiftDist("Jaccard index", p, numthreads), w(1) {}

  void setDist(_shiftDist &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genBitset(p.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genBitset(p.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftCcp : public shiftDist {
 public:
  double w;
 shiftCcp(const Mapfile &p, int numthreads): shiftDist("Cross correlation", p, numthreads), w(1) {}
  
  void setDist(_shiftDist &chr, const vector<char> &fwd, const vector<char> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genVector(p.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genVector(p.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftHamming : public shiftDist {
 public:
  double w;
 shiftHamming(const Mapfile &p, int numthreads): shiftDist("Hamming distance", p, numthreads), w(-1) {}

  void setDist(_shiftDist &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genBitset(p.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genBitset(p.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);
    
    setDist(chr[i], fwd, rev);
  }
};

#endif /* _PW_SHIFTPROFILE_P_H_ */
