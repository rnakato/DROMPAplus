/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_SHIFTPROFILE_P_H_
#define _PW_SHIFTPROFILE_P_H_

#include "pw_gv.h"

std::vector<char> genVector(const strandData &seq, int start, int end);
boost::dynamic_bitset<> genBitset(const strandData &seq, int, int);
void addmp(std::map<int, double> &, const std::map<int, double> &, double w=1);

class ReadShiftProfile {
  int lenF3;
 public:
  std::map<int, double> mp;
  std::map<int, double> nc;
  int start;
  int end;
  int width;
  long nread;
  long len;
  double r;
  double bk;
  double nsc;
  int nsci;
  double rchr;
  
 ReadShiftProfile(const int len, int s=0, int e=0, long n=0, long l=0): lenF3(len), start(s), end(e), width(e-s), nread(n), len(l), r(0), bk(0), nsc(0), nsci(0), rchr(1) {}
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
  void setflen(double w);

  void outputmp(const std::string filename, std::string name, double weight) {
    setControlRatio();
    setflen(weight);
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
  void add2genome(const ReadShiftProfile &x, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
#ifdef DEBUG
    std::cout << "add2genome.." << std::flush;
#endif
    addmp(genome.mp, x.mp, x.rchr);
    addmp(genome.nc, x.nc, x.rchr);
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

    int flen(170);
    int readlen(p.getlenF3());

    std::vector<char> fragarray(chr[i].end - chr[i].start,0);
    std::vector<char> reparray(chr[i].end - chr[i].start,0);
    for(int j=chr[i].start; j<chr[i].end-flen; ++j) {
      if(fwd[j] && rev[j+flen]) for(int k=0; k<flen; ++k) ++fragarray[j - chr[i].start +k];
      if(fwd[j] && rev[j+readlen])  for(int k=0; k<readlen; ++k)  ++reparray[j - chr[i].start +k];
    }

    for(int j=chr[i].start; j<chr[i].end; ++j) {
      if(fragarray[j]<=1) {
	fwd.reset(j);
	rev.reset(j);
      }
    }

    std::vector<range> vrep;
    for(int j=chr[i].start; j<chr[i].end; ++j) {
      if(reparray[j]>=10) j = getRepeatRegion(vrep, j, reparray, chr[i].start, chr[i].end);
    }
    for(auto x:vrep) {
      for(int j=x.start; j<x.end; ++j) {
	fwd.reset(j);
	rev.reset(j);
      }
    }
    
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
