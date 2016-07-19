
/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_SHIFTPROFILE_P_H_
#define _PW_SHIFTPROFILE_P_H_

#include "pw_gv.h"

#define MP_FROM 500
#define MP_TO   1500
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
    int threwidth(5);
    nsc = mp[MP_TO-1]*w;
    for(int i=MP_TO-1-threwidth; i > lenF3*1.3; --i) {
      int on(1);
      for(int j=1; j<=threwidth; ++j) {
	if (mp[i] < mp[i+j] || mp[i] < mp[i-j]) on=0;
      }
      if(on && nsc < mp[i]*r*w) {
	nsc  = mp[i]*r*w;
	nsci = i;
      }
    }
  }

  void outputmp(const string filename, string name, double weight) {
    setControlRatio();
    setflen(weight);
    double sum(getmpsum());
    double rRPKM = (NUM_10M/(double)nread) / (NUM_100M/(double)len);
    double be(bk * rRPKM);
    double const_bu(1/39.0);  // N/(4*L-N), N=10M, L=100M

    ofstream out(filename);
    out << "NSC\t" << nsc*weight << endl;
    out << "RLSC\t"<< (mp[lenF3]*r*weight) << endl;
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
	genome.len   += x.getlenmpbl();
      }
    }
    for(auto x:p.chr) {
      _shiftDist v(p, 0, x.getlen(), x.bothnread_nonred(), x.getlenmpbl());
      //_shiftDist v(p, 180000000, 200000000, x.bothnread_nonred(), x.getlenmpbl());
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

int getRepeatRegion(vector<range> &vrep, int j, vector<int>, int, int);

class shiftJacVec : public shiftDist {
 public:
  double w;
 shiftJacVec(const Mapfile &p, int numthreads): shiftDist("Jaccard index", p, numthreads), w(1) {}

  void setDist(_shiftDist &chr, const vector<char> &fwd, const vector<char> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genVector(p.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genVector(p.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    /*    int segmentsize(1000);
    int nseg = (chr[i].end - chr[i].start)/segmentsize + 1;

    for(int j=0; j<nseg-1; j+=2) {
      for(int k=0; k<segmentsize; ++k) {
	int posi(chr[i].start + j*segmentsize+k);
	fwd[posi] = rev[posi] = 0;
      }
      }*/
    /*    vector<int> fragarray(p.chr[i].getlen(),0);
    vector<int> reparray(p.chr[i].getlen(),0);

    vector<range> vrep;
    int thre=10;
    
    for(int j=chr[i].start; j<chr[i].end; ++j) {
      if(fwd[j] && rev[j+150])          for(int k=0; k<150; ++k)          ++fragarray[j+k];
      if(fwd[j] && rev[j+p.dist.lenF3]) for(int k=0; k<p.dist.lenF3; ++k) ++reparray[j+k];
    }
    for(int j=chr[i].start; j<chr[i].end; ++j) {
      //if( reparray[j]>2 || fragarray[j]>2) cout << j << "\t" << fragarray[j] << "\t" << reparray[j] << endl;
      if(reparray[j]>=thre) j = getRepeatRegion(vrep, j, reparray, chr[i].start, chr[i].end);
    }

    for(auto x:vrep) {
      cout << x.start<< "-" << x.end << endl;
      for(int j=x.start; j<x.end; ++j) fwd[j] = rev[j] = 0; //cout << j<<"\t"<< (int)fwd[j] << "\t" << (int)rev[j] << endl;  //
      }*/

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
