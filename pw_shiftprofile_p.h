/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_SHIFTPROFILE_P_H_
#define _PW_SHIFTPROFILE_P_H_

#include "pw_gv.h"

namespace {
  const int mp_from(500);
  const int mp_to(1500);
  const int ng_from(4000);
  const int ng_to(5000);
  const int ng_step(100);
}

std::vector<char> genVector(const strandData &seq, int start, int end);
boost::dynamic_bitset<> genBitset(const strandData &seq, int, int);
void addmp(std::map<int, double> &, const std::map<int, double> &, double w=1);

class FragmentVariability {
  int fraglen;
  long numOfFragment;
  long numOfCoveredBase;
  long nread;
  std::vector<int> FragOverlapDist;
  static const int SizeOfFragOverlapDist=10000;
  double SumOfFragOverlapDist;
  
 public:
 FragmentVariability(): fraglen(0), numOfFragment(0), numOfCoveredBase(0), nread(0), FragOverlapDist(SizeOfFragOverlapDist,0), SumOfFragOverlapDist(0) {}
  void setVariability(const int l, const int start, const int end ,const int width,
		      const boost::dynamic_bitset<> &fwd,
		      const boost::dynamic_bitset<> &rev) {
    fraglen = l;
    std::vector<short> array(width, 0);
    for(int i=start; i<end - fraglen; ++i) {
      if(fwd[i] && rev[i+fraglen]) {
	++numOfFragment;
	for(int j=0; j<fraglen; ++j) ++array[i - start +j];
      }
    }
    for(int i=start; i<end - fraglen; ++i) if(array[i]) ++numOfCoveredBase;
    for(int i=start; i<end; ++i) ++FragOverlapDist[array[i]];
    SumOfFragOverlapDist = accumulate(FragOverlapDist.begin(), FragOverlapDist.end(), 0);
  }
  double getFragOverlapDist(const int i) const {
    if(i<0 || i>= SizeOfFragOverlapDist) {
      std::cerr << "error: invalid num " << i << "for getFragOverlapDist. max: " << SizeOfFragOverlapDist << std::endl;
      return -1;
    }
    else return FragOverlapDist[i]/SumOfFragOverlapDist;
  }
  double getUniqueness() const {
    return numOfCoveredBase / (double)(numOfFragment * fraglen);
  }
  void add2genome(const FragmentVariability &x) {
    fraglen           = x.fraglen;
    numOfFragment    += x.numOfFragment;
    numOfCoveredBase += x.numOfCoveredBase;
    nread            += x.nread;
    for(uint i=0; i<FragOverlapDist.size(); ++i) FragOverlapDist[i] += x.FragOverlapDist[i];
    SumOfFragOverlapDist = accumulate(FragOverlapDist.begin(), FragOverlapDist.end(), 0);
  }
};

class ReadShiftProfile {
  int lenF3;
  double nsc;
  int nsci;
  double r;
  double bk;
  
 public:
  std::map<int, double> mp;
  std::map<int, double> nc;
  int start;
  int end;
  int width;
  long nread;

  long len;
  double rchr;

  FragmentVariability fvfrag;
  FragmentVariability fvrep;
  FragmentVariability fvback;

 ReadShiftProfile(const int len, int s=0, int e=0, long n=0, long l=0): lenF3(len), nsc(0), nsci(0), r(0), bk(0), start(s), end(e), width(e-s), nread(n), len(l), rchr(1) {}
  void setmp(int i, double val, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    mp[i] = val;
  }
  void setFragmentVariability4Frag(const int l, const boost::dynamic_bitset<> &fwd, const boost::dynamic_bitset<> &rev) {
    fvfrag.setVariability(l, start, end, width, fwd, rev);
  }
  void setFragmentVariability4Rep(const int l, const boost::dynamic_bitset<> &fwd, const boost::dynamic_bitset<> &rev) {
    fvrep.setVariability(l, start, end, width, fwd, rev);
  }
  void setFragmentVariability4Back(const int l, const boost::dynamic_bitset<> &fwd, const boost::dynamic_bitset<> &rev) {
    fvback.setVariability(l, start, end, width, fwd, rev);
  }
  
  int getnsci() const { return nsci; }
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

  void outputmp(const std::string filename, std::string name, double weight) {
    double sum(getmpsum());
    double rRPKM = (NUM_10M/static_cast<double>(nread)) / (NUM_100M/static_cast<double>(len));
    double be(bk * rRPKM);
    double const_bu(1/39.0);  // N/(4*L-N), N=10M, L=100M

    std::ofstream out(filename);
    out << "NSC\t" << nsc*weight << std::endl;
    out << "RLSC\t"<< (mp[lenF3]*r*weight) << std::endl;
    out << "Estimated fragment length\t" << nsci << std::endl;
    out << "Background enrichment\t" << be << std::endl;
    out << "Background uniformity\t" << const_bu / be << std::endl;
    //    out << "nread_peak\t" << nread_peak << "\t" << nread_peak/(double)nread*100 << std::endl;
    // out << "nread_back\t" << nread_back << "\t" << nread_back/(double)nread*100 << std::endl;
    // out << "nread_rep\t"  << nread_rep  << "\t" << nread_rep/(double)nread*100  << std::endl;
    out << "Fraglen uniqueness\t" << fvfrag.getUniqueness() << std::endl;
    out << "Readlen uniqueness\t" << fvrep.getUniqueness() << std::endl;
    out << "Background uniqueness\t" << fvback.getUniqueness() << std::endl;
    out << "Normalized Fraglen uniqueness\t" << fvfrag.getUniqueness()/fvback.getUniqueness() << std::endl;

    out << "Strand shift\t" << name << "\tprop\tper 10M reads\tper control" << std::endl;
    for(auto itr = mp.begin(); itr != mp.end(); ++itr) 
      out << itr->first            << "\t"
	  << itr->second           << "\t"
	  << (itr->second/sum)     << "\t"
	  << (itr->second * rRPKM) << "\t"
	  << (itr->second * r)     << std::endl;
  }
};

class ReadShiftProfileAll {
 protected:
  std::vector<range> seprange;
  
 public:
  std::string name;
  ReadShiftProfile genome;
  std::vector<ReadShiftProfile> chr;
  
  void defSepRange(int numthreads);
 ReadShiftProfileAll(std::string n, const Mapfile &p, int numthreads): name(n), genome(p.getlenF3()) {
    for(auto x:p.chr) {
      if(x.isautosome()) {
	genome.nread += x.bothnread_nonred();
	genome.len   += x.getlenmpbl();
      }
    }
    for(auto x:p.chr) {
      ReadShiftProfile v(p.getlenF3(), 0, x.getlen(), x.bothnread_nonred(), x.getlenmpbl());
      //ReadShiftProfile v(p, 0, 120000000, x.bothnread_nonred(), x.getlenmpbl());
      v.rchr = v.nread/static_cast<double>(genome.nread);
      chr.push_back(v);
    }
    // seprange
    defSepRange(numthreads);
  }
  void addmp2genome(const ReadShiftProfile &x) {
    addmp(genome.mp, x.mp, x.rchr);
    addmp(genome.nc, x.nc, x.rchr);
  }
  void addnread2genome(const ReadShiftProfile &x, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    genome.fvfrag.add2genome(x.fvfrag);
    genome.fvrep.add2genome(x.fvrep);
    genome.fvback.add2genome(x.fvback);
  }
};

int getRepeatRegion(std::vector<range> &vrep, int j, std::vector<char>, int, int);

class shiftJacVec : public ReadShiftProfileAll {
 public:
  double w;
 shiftJacVec(const Mapfile &p, int numthreads): ReadShiftProfileAll("Jaccard index", p, numthreads), w(1) {}

  void setDist(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genVector(p.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genVector(p.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    /*    for(int j=chr[i].start; j<chr[i].end; ++j) {
      if(fwd[j] && rev[j+170]) std::cout << j << "\t" << fwd[j] << "\t" << rev[j+170] << std::endl;
      }*/
    
    /*    std::vector<int> fragarray(p.chr[i].getlen(),0);
    std::vector<int> reparray(p.chr[i].getlen(),0);

    //    std::vector<range> vrep;
    int thre=3;
    
    for(int j=chr[i].start; j<chr[i].end-170; ++j) {
      if(fwd[j] && rev[j+170]) for(int k=0; k<170; ++k) ++fragarray[j+k];
      if(fwd[j] && rev[j+50])  for(int k=0; k<50; ++k)  ++reparray[j+k];
    }

    std::vector<int> dfrag(170,0);
    std::vector<int> drep(170,0);
    for(int j=chr[i].start; j<chr[i].end; ++j) {
      ++dfrag[fragarray[j]];
      ++drep[reparray[j]];
    }
    int ndfragon(0);
    int ndrepon(0);
    for(int j=chr[i].start; j<chr[i].end; ++j) {
      //      if(fragarray[j]>=thre) j = getRepeatRegion(vrep, j, reparray, chr[i].start, chr[i].end);
      if(fragarray[j]) ++ndfragon;
      if(reparray[j]) ++ndrepon;
    }
    for(int j=0; j<170; ++j) std::cout << j << "\t" << dfrag[j] << "\t" << drep[j] << std::endl;
    

    std::cout <<"covered num: " << ndfragon << "\t" << ndrepon << std::endl;
    
    for(int j=chr[i].start; j<chr[i].end; ++j) {
      if(fragarray[j]<=1) fwd[j] = rev[j] = 0;
      }*/
    /*    for(auto x:vrep) {
      std::cout << x.start<< "-" << x.end << std::endl;
      for(int j=x.start; j<x.end; ++j) fwd[j] = rev[j] = 0;
      }*/

    setDist(chr[i], fwd, rev);  
  }
};

class shiftJacBit : public ReadShiftProfileAll {
 public:
  double w;
 shiftJacBit(const Mapfile &p, int numthreads): ReadShiftProfileAll("Jaccard index", p, numthreads), w(1) {}

  void setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genBitset(p.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genBitset(p.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftCcp : public ReadShiftProfileAll {
 public:
  double w;
 shiftCcp(const Mapfile &p, int numthreads): ReadShiftProfileAll("Cross correlation", p, numthreads), w(1) {}
  
  void setDist(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genVector(p.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genVector(p.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftHamming : public ReadShiftProfileAll {
 public:
  double w;
 shiftHamming(const Mapfile &p, int numthreads): ReadShiftProfileAll("Hamming distance", p, numthreads), w(-1) {}

  void setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genBitset(p.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genBitset(p.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);
    
    setDist(chr[i], fwd, rev);
  }
};

#endif /* _PW_SHIFTPROFILE_P_H_ */
