/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_SHIFTPROFILE_P_H_
#define _PW_SHIFTPROFILE_P_H_

#include "pw_gv.h"
#include <boost/dynamic_bitset.hpp>

namespace {
  const int mp_from(500);
  const int mp_to(1500);
  const int ng_from(4000);
  const int ng_to(5000);
  const int ng_step(100);
  //  const int SizeOfFragOverlapDist(10000);
  const int sizeOfvDistOfDistaneOfFrag = 2000;
}

std::vector<char> genVector(const strandData &seq, int start, int end);
boost::dynamic_bitset<> genBitset(const strandData &seq, int, int);
void addmp(std::map<int, double> &, const std::map<int, double> &, double w);

class FragmentVariability {
  //  std::vector<int> FragOverlapDist;
  double sumOfvDistOfDistaneOfFrag;
  static const int length4Accu = 500;
  std::vector<int> vDistOfDistaneOfFrag;
  std::vector<double> vAccuOfDistaneOfFrag;

 public:
 FragmentVariability(): sumOfvDistOfDistaneOfFrag(0),
    vDistOfDistaneOfFrag(sizeOfvDistOfDistaneOfFrag, 0),
    vAccuOfDistaneOfFrag(sizeOfvDistOfDistaneOfFrag, 0) {} //, FragOverlapDist(SizeOfFragOverlapDist,0)

  void setVariability(const int fraglen, const int start, const int end, const int width,
		      const std::vector<char> &fwd, const std::vector<char> &rev) {
    int s(start);
    if(start + fraglen < 0) s = start - fraglen;
    int e(end - std::max(fraglen, length4Accu));
    int last(s);
    for(int i=s; i<e; ++i) {
      if(fwd[i] && rev[i+fraglen]) {
	if(RANGE(i-last, 0, sizeOfvDistOfDistaneOfFrag-1)) ++vDistOfDistaneOfFrag[i-last];
	last = i;
      }
    }

    sumOfvDistOfDistaneOfFrag = accumulate(vDistOfDistaneOfFrag.begin(), vDistOfDistaneOfFrag.end(), 0);

    double a(0);
    for(int i=0; i<sizeOfvDistOfDistaneOfFrag; ++i) {
      a += getDistOfDistanceOfFragment(i);
      vAccuOfDistaneOfFrag[i] = a;
    }
  }

  //  long getseqsize() const { return accumulate(FragOverlapDist.begin(), FragOverlapDist.end(), 0); }
  /*  double getPropOfOverlapDist(const int i) const {
    if(i<0 || i>= SizeOfFragOverlapDist) {
      std::cerr << "error: invalid num " << i << "for getPropOfOverlapDist. max: " << SizeOfFragOverlapDist << std::endl;
      return -1;
    }
    else {
      return FragOverlapDist[i]/static_cast<double>(accumulate(FragOverlapDist.begin(), FragOverlapDist.end(), 0));
    }
    }*/
  double getDistOfDistanceOfFragment(const int i) const {
    if(i<0 || i>= sizeOfvDistOfDistaneOfFrag) {
      std::cerr << "error: invalid num " << i << "for getDistOfDistanceOfFragment max: " << sizeOfvDistOfDistaneOfFrag << std::endl;
      return -1;
    }
    else {
      return vDistOfDistaneOfFrag[i] / sumOfvDistOfDistaneOfFrag;
    }
  }
  double getAccuOfDistanceOfFragment(const int i) const {
    if(i<0 || i>= sizeOfvDistOfDistaneOfFrag) {
      std::cerr << "error: invalid num " << i << "for getAccuOfDistanceOfFragment max: " << sizeOfvDistOfDistaneOfFrag << std::endl;
      return -1;
    }
    else {
      return vAccuOfDistaneOfFrag[i];
    }
  }
  void add2genome(const FragmentVariability &x) {
    //    for(uint i=0; i<FragOverlapDist.size(); ++i) FragOverlapDist[i] += x.FragOverlapDist[i];
    for(uint i=0; i<vDistOfDistaneOfFrag.size(); ++i) vDistOfDistaneOfFrag[i] += x.vDistOfDistaneOfFrag[i];
  }
};

class ReadShiftProfile {
  int lenF3;
  double nsc;
  int nsci;
  double r;
  double bk;
  double w;

 protected:
  long len;
  long nread;
  
 public:
  std::map<int, double> mp;
  std::map<int, double> nc;
  int start;
  int end;
  int width;

  double rchr;

 ReadShiftProfile(const int len, const double wref, int s=0, int e=0, long n=0, long l=0): lenF3(len), nsc(0), nsci(0), r(0), bk(0), w(wref), len(l), nread(n), start(s), end(e), width(e-s), rchr(1) {}
  virtual ~ReadShiftProfile() {}
  void setmp(int i, double val, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    mp[i] = val;
  }

  void setrchr(const long n) { rchr = n ? nread/static_cast<double>(n): 0; }
  int getnsci() const { return nsci; }
  double getmpsum() const {
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
    int threwidth(5);
    setControlRatio();
    nsc = mp[mp_to-1]*w;
    for(int i=mp_to-1-threwidth; i > lenF3*1.3; --i) {
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

  void print2file(const std::string filename, const std::string name) const {
    if(!nread) {
      std::cerr << filename << ": no read" << std::endl;
    }
    double sum(getmpsum());
    double rRPKM = (NUM_10M/static_cast<double>(nread)) / (NUM_100M/static_cast<double>(len));
    double be(bk * rRPKM);
    double const_bu(1/39.0);  // N/(4*L-N), N=10M, L=100M

    std::ofstream out(filename);
    out << "NSC\t" << nsc*w << std::endl;
    out << "RLSC\t"<< (mp.at(lenF3)*r*w) << std::endl;
    out << "Estimated fragment length\t" << nsci << std::endl;
    out << "Background enrichment\t" << be << std::endl;
    out << "Background uniformity\t" << const_bu / be << std::endl;

    out << "Strand shift\t" << name << "\tprop\tper 10M reads\tper control" << std::endl;
    for(auto itr = mp.begin(); itr != mp.end(); ++itr) 
      out << itr->first            << "\t"
	  << itr->second           << "\t"
	  << (itr->second/sum)     << "\t"
	  << (itr->second * rRPKM) << "\t"
	  << (itr->second * r)     << std::endl;
  }
};

class ReadShiftProfileGenome: public ReadShiftProfile {
  std::string name;

 protected:
  std::vector<range> seprange;
  
 public:
  std::vector<ReadShiftProfile> chr;
  
 ReadShiftProfileGenome(std::string n, const Mapfile &p, const int numthreads, double wref): name(n), ReadShiftProfile(p.getlenF3(), wref) {
    for(auto x:p.genome.chr) {
      if(x.isautosome()) {
	nread += x.bothnread_nonred();
	len   += x.getlenmpbl();
      }
    }
    for(auto x:p.genome.chr) {
      ReadShiftProfile v(p.getlenF3(), wref, 0, x.getlen(), x.bothnread_nonred(), x.getlenmpbl());
      //ReadShiftProfile v(p, 0, 120000000, x.bothnread_nonred(), x.getlenmpbl());
      v.setrchr(nread);
      chr.push_back(v);
    }
    // seprange
    defSepRange(numthreads);
  }
  virtual ~ReadShiftProfileGenome(){}
  void defSepRange(const int numthreads) {
    int length(mp_to+mp_from);
    int sepsize = length/numthreads +1;
    for(int i=0; i<numthreads; ++i) {
      int s = i*sepsize;
      int e = (i+1)*sepsize;
      if(i==numthreads-1) e = length;
      seprange.push_back(range(s - mp_from, e - mp_from));
    }
  }
  void addmp2genome(const int i) {
    addmp(mp, chr[i].mp, chr[i].rchr);
    addmp(nc, chr[i].nc, chr[i].rchr);
  }
  void outputmpGenome(const std::string &filename) const {
    print2file(filename, name);
  }
  void outputmpChr(const std::string &filename, const int i) const {
    chr[i].print2file(filename, name);
  }
  void printStartMessage() const {
    std::cout << "Making " << name << " profile..." << std::flush;
  }
  
};

int getRepeatRegion(std::vector<range> &vrep, int j, std::vector<char>, int, int);

class shiftJacVec : public ReadShiftProfileGenome {
 public:
 shiftJacVec(const Mapfile &p, int numthreads): ReadShiftProfileGenome("Jaccard index", p, numthreads, 1) {}

  void setDist(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genVector(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genVector(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);  
  }
};

class shiftJacBit : public ReadShiftProfileGenome {
 public:
 shiftJacBit(const Mapfile &p, int numthreads): ReadShiftProfileGenome("Jaccard index", p, numthreads, 1) {}

  void setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genBitset(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genBitset(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftCcp : public ReadShiftProfileGenome {
 public:
 shiftCcp(const Mapfile &p, int numthreads): ReadShiftProfileGenome("Cross correlation", p, numthreads, 1) {}
  
  void setDist(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genVector(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genVector(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftHamming : public ReadShiftProfileGenome {
 public:
 shiftHamming(const Mapfile &p, int numthreads): ReadShiftProfileGenome("Hamming distance", p, numthreads, -1) {}

  void setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genBitset(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genBitset(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);
    
    setDist(chr[i], fwd, rev);
  }
};

class shiftFragVar : public ReadShiftProfileGenome {
 public:
 shiftFragVar(const Mapfile &p, int numthreads): ReadShiftProfileGenome("Fragment Variability", p, numthreads, 1) {}

  void setDist(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genVector(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genVector(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

#endif /* _PW_SHIFTPROFILE_P_H_ */
