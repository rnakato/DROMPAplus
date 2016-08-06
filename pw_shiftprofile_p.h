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
  const int sizeOfvDistOfDistaneOfFrag = 5000;

  const std::vector<int> v4mpfv{50, 100, 150, 500, 1000, 2000, 3000, 10000, 100000, 1000000};
}

std::vector<char> genVector(const strandData &seq, int start, int end);
std::vector<char> genVector4FixedReadsNum(const strandData &seq, int start, int end, const double r4cmp);
boost::dynamic_bitset<> genBitset(const strandData &seq, int, int);
void addmp(std::map<int, double> &, const std::map<int, double> &, double w);

class FragmentVariability {
  double sumOfvDistOfDistaneOfFrag;
  std::vector<int> vDistOfDistaneOfFrag;

 public:
 FragmentVariability(): sumOfvDistOfDistaneOfFrag(0),
    vDistOfDistaneOfFrag(sizeOfvDistOfDistaneOfFrag, 0) {}

  void setVariability(const int fraglen, const int start, const int end,
		      const std::vector<char> &fwd, const std::vector<char> &rev) {
    int s(start);
    if(start + fraglen < 0) s = start - fraglen;
    int e(end - fraglen);
    int last(s);
    for(int i=s; i<e; ++i) {
      if(fwd[i] && rev[i+fraglen]) {
	if(RANGE(i-last, 0, sizeOfvDistOfDistaneOfFrag-1)) ++vDistOfDistaneOfFrag[i-last];
	else ++vDistOfDistaneOfFrag[sizeOfvDistOfDistaneOfFrag-1];
	last = i;
      }
    }

    sumOfvDistOfDistaneOfFrag = accumulate(vDistOfDistaneOfFrag.begin(), vDistOfDistaneOfFrag.end(), 0);
  }

  double getDistOfDistanceOfFragment(const int i) const {
    if(i<0 || i>= sizeOfvDistOfDistaneOfFrag) {
      std::cerr << "error: invalid num " << i << "for getDistOfDistanceOfFragment max: " << sizeOfvDistOfDistaneOfFrag << std::endl;
      return -1;
    }
    else {
      return sumOfvDistOfDistaneOfFrag ? vDistOfDistaneOfFrag[i] / sumOfvDistOfDistaneOfFrag : 0;
    }
  }
  double getAccuOfDistanceOfFragment(const int i) const {
    if(i<0 || i>= sizeOfvDistOfDistaneOfFrag) {
      std::cerr << "error: invalid num " << i << "for getAccuOfDistanceOfFragment max: " << sizeOfvDistOfDistaneOfFrag << std::endl;
      return -1;
    }
    else {
      double vAccuOfDistaneOfFrag(0);
      for(int j=0; j<=i; ++j) {
	vAccuOfDistaneOfFrag += getDistOfDistanceOfFragment(j);
      }
      return vAccuOfDistaneOfFrag;
    }
  }
  void add2genome(const FragmentVariability &x, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    for(uint i=0; i<vDistOfDistaneOfFrag.size(); ++i) vDistOfDistaneOfFrag[i] += x.vDistOfDistaneOfFrag[i];
    sumOfvDistOfDistaneOfFrag += x.sumOfvDistOfDistaneOfFrag;
  }
};

class ReadShiftProfile {
  int lenF3;
  double nsc;
  int nsci;
  double r;
  double bk;
  double dirOfProfileSummit;

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

 ReadShiftProfile(const int len, const double wref, int s=0, int e=0, long n=0, long l=0):
  lenF3(len), nsc(0), nsci(0), r(0), bk(0), dirOfProfileSummit(wref), len(l), nread(n), start(s), end(e), width(e-s), rchr(1) {}
  virtual ~ReadShiftProfile() {}
  void setmp(const int i, const double val, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    mp[i] = val;
  }

  void setrchr(const long n) { rchr = n ? nread/static_cast<double>(n): 0; }
  int getlenF3() const { return lenF3; }
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
    nsc = mp[mp_to-1] * dirOfProfileSummit;
    for(int i=mp_to-1-threwidth; i > lenF3*1.3; --i) {
      int on(1);
      for(int j=1; j<=threwidth; ++j) {
	if (mp[i] < mp[i+j] || mp[i] < mp[i-j]) on=0;
      }
      if(on && nsc < mp[i] *r *dirOfProfileSummit) {
	nsc  = mp[i] *r *dirOfProfileSummit;
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
    out << "NSC\t" << nsc * dirOfProfileSummit << std::endl;
    out << "RLSC\t"<< (mp.at(lenF3) *r *dirOfProfileSummit) << std::endl;
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
  void print2file4fvp(const std::string filename, const std::string name, const int flen, const bool lackOfReads) const {
    if(!nread) {
      std::cerr << filename << ": no read" << std::endl;
    }
    std::ofstream out(filename);

    //    std::cout << flen <<"," << lenF3 <<std::endl;

    std::string str("");
    if(lackOfReads) str = " (read number is insufficient)";
    out << "Fragment score" << str << "\t" << mp.at(flen) << std::endl;
    out << "Read score" << str << "\t" << mp.at(lenF3) << std::endl;
    out << "Strand shift\t" << name << std::endl;
    for(auto itr = mp.begin(); itr != mp.end(); ++itr) 
      out << itr->first << "\t" << itr->second << std::endl;
  }
};

class ReadShiftProfileGenome: public ReadShiftProfile {
 protected:
  std::string name;
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
      v.setrchr(nread);
      chr.push_back(v);
    }
    // seprange
    defSepRange(numthreads);
  }
  virtual ~ReadShiftProfileGenome(){}
  long getnread() const { return nread; }
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
  
  void printmpfv(const std::string &){};
};

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
  std::map<int, FragmentVariability> mpfv;
  int flen;
  bool lackOfReads;
  bool fvpfull;
 public:
 shiftFragVar(const Mapfile &p, const int numthreads, const int fl, const bool b):
  ReadShiftProfileGenome("Fragment Variability", p, numthreads, 1), flen(fl), lackOfReads(false), fvpfull(b) {}

  void setDist(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev);
  void execchr(const Mapfile &p, const int i, const double r4cmp) {
    auto fwd = genVector4FixedReadsNum(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end, r4cmp);
    auto rev = genVector4FixedReadsNum(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end, r4cmp);

    //    scanRepeatRegion(fwd, rev);
    
    setDist(chr[i], fwd, rev);
  }

  void lackOfReads_on() { lackOfReads=true; }
  void printmpfv(const std::string &filename) const {
    std::ofstream out(filename);

    /*    out << "Accumulated: " << std::endl;
    for(auto x: v4mpfv) {
      if(fvpfull && x > mp_to) continue;
      double sum(0);
      for(size_t k=0; k<sizeOfvDistOfDistaneOfFrag; ++k) sum += mpfv.at(x).getAccuOfDistanceOfFragment(k);
      out << "len" << x << "\t" << sum << std::endl;
      }*/
    
    for(auto x: v4mpfv) {
      if(fvpfull && x > mp_to) continue;
      out << "\tlen" << x;
    }
    for(auto x: v4mpfv) {
      if(fvpfull && x > mp_to) continue;
      out << "\tlen" << x;
    }
    out << std::endl;
    for(size_t k=0; k<sizeOfvDistOfDistaneOfFrag; ++k) {
      out << k << "\t";

      for(auto x: v4mpfv) {
	if(fvpfull && x > mp_to) continue;
	out << mpfv.at(x).getAccuOfDistanceOfFragment(k) << "\t";
      }
      for(auto x: v4mpfv) {
	if(fvpfull && x > mp_to) continue;
	out <<  mpfv.at(x).getDistOfDistanceOfFragment(k) << "\t";
      }
      out << std::endl;
    }
  }
  
  void outputmpGenome(const std::string &filename) const {
    print2file4fvp(filename, name, flen, lackOfReads);
  }
  void outputmpChr(const std::string &filename, const int i) const {  
    chr[i].print2file4fvp(filename, name, flen, lackOfReads);
  }
};

#endif /* _PW_SHIFTPROFILE_P_H_ */
