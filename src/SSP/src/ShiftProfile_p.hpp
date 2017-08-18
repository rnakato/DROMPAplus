/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _SSP_SHIFTPROFILE_P_H_
#define _SSP_SHIFTPROFILE_P_H_

#include <boost/dynamic_bitset.hpp>
#include "ssp_gv.hpp"
#include "../common/alglib/alglib.h"

namespace {
  const int32_t mp_from(500);
  const int32_t mp_to(1500);
}

std::vector<int8_t> genVector(const std::vector<Read> &vReadref, const int32_t start, const int32_t end);
boost::dynamic_bitset<> genBitset(const std::vector<Read> &vReadref, const int32_t, const int32_t);
void addmp(std::map<int32_t, double> &, const std::map<int32_t, double> &, double w);

double getmpmean(const std::map<int32_t, double> mp, int32_t s, int32_t e) {
  double m(0);
  int32_t n(0);
  if(s > e) {
    std::cerr << "Error: s " << s << "> e " << e <<"for getmpmean" << std::endl;
    return -1;
  }
  for(int32_t i=s; i <= e; ++i) {
    m += mp.at(i);
    ++n;
  }
  return m/n;
}

class ReadShiftProfile {
  int32_t lenF3;
  double r;
  double bk;
  int32_t bk_from;

 protected:
  double nsc;
  double rsc;
  double rlsc;
  int32_t nsci;
  uint64_t len;
  uint64_t nread;
  uint32_t num4ssp;
  double backgroundUniformity;
  
 public:
  std::map<int32_t, double> mp;
  std::map<int32_t, double> nc;
  int32_t start;
  int32_t end;
  int32_t width;

  double rchr;

 ReadShiftProfile(const int32_t lenf3, const int32_t b, const int32_t n4s, int32_t s=0, int32_t e=0, int64_t n=0, int64_t l=0):
  lenF3(lenf3), r(0), bk(0), bk_from(b), nsc(0), rsc(0), rlsc(0), nsci(0), len(l), nread(n), num4ssp(n4s), backgroundUniformity(0), start(s), end(e), width(e-s), rchr(1) {}
  virtual ~ReadShiftProfile() {}

  void setmp(const int32_t i, const double val, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    mp[i] = val;
  }

  double getbackgroundUniformity() const { return backgroundUniformity; }
  void setrchr(const uint64_t n) { rchr = n ? getratio(nread, n): 0; }
  int32_t getlenF3() const { return lenF3; }
  double getnsc()    const { return nsc; }
  double getrlsc()    const { return rlsc; }
  double getrsc()    const { return rsc; }
  int32_t getnsci()  const { return nsci; }
  uint64_t getnread() const { return nread; }
  uint64_t getlen()  const { return len; }
  double getmpsum()  const {
    double sum(0);
    for(auto pair: mp) sum += pair.second;
    return sum;
  }
  void setControlRatio() {
    int32_t n(0);
    for(auto pair: nc) {
      if(pair.first >= bk_from) {
	bk += pair.second;
	++n;
      }
    }
    bk /= n;
    r = 1/bk;
  }

  void setflen(const std::string &name) {
    int32_t threwidth(5);
    int32_t leftend(lenF3*1.2);
    if(leftend>150) leftend=150;

    setControlRatio();
    if(name == "Hamming distance") nsc = mp[mp_to-1];
    else nsc = mp[mp_to-1]*r;

    std::map<int32_t, double> mpsmooth;
    for(int32_t i=mp_to-1-2; i > leftend-threwidth; --i) {
      mpsmooth[i] = getmpmean(mp, i-2, i+2);
    }

    for(int32_t i=mp_to-1-threwidth-2; i > leftend; --i) {
      if(name == "Hamming distance") {
	if(mpsmooth.at(i) > mpsmooth.at(i+threwidth) || mpsmooth.at(i) > mpsmooth.at(i-threwidth)) continue;
	if(nsc > mp.at(i)) {
	  nsc  = mp.at(i);
	  nsci = i;
	}
      } else {
	if(mpsmooth.at(i) < mpsmooth.at(i+threwidth) || mpsmooth.at(i) < mpsmooth.at(i-threwidth)) continue;
	double s(mp.at(i)*r);
	if(nsc < s) {
	  nsc  = s;
	  rsc  = (mp.at(i) - bk)/(mp.at(lenF3) - bk);
	  nsci = i;
	}
      }
    }
  }

  void print2file(const std::string &filename, const std::string &name) {
    if(!nread) {
      std::cerr << filename << ": no read" << std::endl;
    }
    double sum(getmpsum());
    double rRPKM = getratio(num4ssp, nread) / getratio(NUM_100M, len);

    double be(bk * rRPKM);
    double const_bu = getratio(num4ssp, (4*NUM_100M - num4ssp));  // 1/39 N/(4*L-N), N=10M, L=100M
    //    std::cout << "####### " << num4ssp << "\t  " << const_bu << "\t" << rRPKM << "\t" << bk << std::endl;
    //    std::cout << "####### " << nread << "\t  " << NUM_100M << "\t" << len << std::endl;
    rlsc = mp.at(lenF3) *r;
    backgroundUniformity = const_bu / be;

    std::ofstream out(filename);
    out << "NSC\t" << nsc << std::endl;
    out << "RSC\t" << rsc << std::endl;
    out << "RLSC\t"<< rlsc << std::endl;
    out << "Estimated fragment length\t" << nsci << std::endl;
    out << "Background enrichment\t" << be << std::endl;
    out << "Background uniformity\t" << backgroundUniformity << std::endl;

    out << "Strand shift\t" << name << "\tProportion\tper " << num4ssp/NUM_1M << "M reads for 100Mbp len\tper control" << std::endl;
    for(auto pair: mp)
      out << pair.first            << "\t"
	  << pair.second           << "\t"
	  << (pair.second/sum)     << "\t"
	  << (pair.second * rRPKM) << "\t"
	  << (pair.second * r)     << std::endl;
    for(auto pair: nc)
      out << pair.first            << "\t"
	  << pair.second           << "\t"
	  << (pair.second/sum)     << "\t"
	  << (pair.second * rRPKM) << "\t"
	  << (pair.second * r)     << std::endl;
  }
};

class ReadShiftProfileGenome: public ReadShiftProfile {
 protected:
  int32_t ng_from;
  int32_t ng_to;
  int32_t ng_step;
  std::vector<range> seprange;
  
 public:
  std::string name;
  std::vector<ReadShiftProfile> chr;
  
 ReadShiftProfileGenome(const std::string n, const SSPstats &sspst, const SeqStatsGenome &genome):
  ReadShiftProfile(genome.dflen.getlenF3(), sspst.getNgFrom(), sspst.getnum4ssp()),
    ng_from(5000),
    ng_to(sspst.getNgTo()),
    ng_step(sspst.getNgStep()),
    name(n)
    {
      for(auto &x: genome.chr) {
	if(x.isautosome()) {
	  nread += x.getnread_nonred(Strand::BOTH);
	  len   += x.getlenmpbl();
	  //	  std::cout<< len << "\t" << x.getlenmpbl() << std::endl;
	}
      }
      for(auto &x: genome.chr) {
	ReadShiftProfile v(genome.dflen.getlenF3(), sspst.getNgFrom(), sspst.getnum4ssp(), 0, x.getlen(), x.getnread_nonred(Strand::BOTH), x.getlenmpbl());
	v.setrchr(nread);
	chr.push_back(v);
      }
      // seprange
      defSepRange(sspst.getnumthreads());
    }
  virtual ~ReadShiftProfileGenome(){}

  void defSepRange(const int32_t numthreads) {
    int32_t length(mp_to + mp_from);
    int32_t sepsize(length/numthreads +1);
    for(int32_t i=0; i<numthreads; ++i) {
      int32_t s = i*sepsize;
      int32_t e = (i+1)*sepsize;
      if(i==numthreads-1) e = length;
      seprange.push_back(range(s - mp_from, e - mp_from));
    }
  }
  void addmp2genome(const int32_t i) {
    addmp(mp, chr[i].mp, chr[i].rchr);
    addmp(nc, chr[i].nc, chr[i].rchr);
  }

  void makeRscript(const std::string &filename, const std::string &prefix) {
    std::string Rscript(prefix + ".R");
    std::ofstream out(Rscript);

    out << "data <- read.csv('" << filename << "', header=TRUE, skip=6, sep='\\t', quote='')" << std::endl;
    out << "output <- '" << prefix << "'" << std::endl;
    out << "pdf(paste(output, '.pdf', sep=''), height=7, width=14)" << std::endl;
    out << "par(mfrow=c(1,2))" << std::endl;
    out << "plot(data[1:" << mp_from+mp_to << ",1], data[1:" << mp_from+mp_to << ",5], type='l', xlab='Strand shift', ylab='Score relative to background', xlim=c(" << -mp_from << "," << mp_to << "), main='" << -mp_from << " bp ~ " << mp_to << " bp', ";
    if(name == "Jaccard index") out << "sub=sprintf('NSC=%g, RSC=%g, RLSC=%g, Bu=%g', " << nsc << "," << rsc << ","  << rlsc << "," << backgroundUniformity << "))" << std::endl;
    else if(name == "Cross correlation") out << "sub=sprintf('NSC=%g, RSC=%g, RLSC=%g', " << nsc << "," << rsc << ","  << rlsc << "))" << std::endl;
    else out << ")" << std::endl;
    out << "abline(v=" << nsci <<",lty=2,col=2)" << std::endl;
    out << "abline(v=" << getlenF3() <<",lty=2,col='blue')" << std::endl;
    out << "legend('bottomright', legend=paste('Estimated fragment length = ', " << nsci << "))" << std::endl;
    out << "plot(data[,1], data[,5], type='l', xlab='Strand shift',ylab='Score relative to background', main='Long distance (log scale)', log='x', xlim=c(1,1000000))" << std::endl;
    out << "dev.off()" << std::endl;

    std::string command = "R --vanilla < " + Rscript + " > " + Rscript + ".log 2>&1";
    
    int32_t return_code = system(command.c_str());
    if(WEXITSTATUS(return_code)) {
      std::cerr << "Warning: command " << command << "return nonzero status." << std::endl;
    }
  }
  
  void outputmpGenome(const std::string &prefix) {
    std::string filename = prefix + ".csv";
    print2file(filename, name);
    makeRscript(filename, prefix);
  }
  void outputmpChr(const std::string &filename, const int32_t i) {
    chr[i].print2file(filename, name);
  }
  void printStartMessage() const {
    std::cout << "Making " << name << " profile..." << std::flush;
  }
};

class shiftJacVec : public ReadShiftProfileGenome {
 public:
 shiftJacVec(const SSPstats &sspst, const SeqStatsGenome &genome):
   ReadShiftProfileGenome("Jaccard index", sspst, genome) {}

  void setDist(ReadShiftProfile &chr, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev);
  void execchr(const SeqStatsGenome &genome, int32_t i) {
    auto fwd = genVector(genome.chr[i].getvReadref(Strand::FWD), chr[i].start, chr[i].end);
    auto rev = genVector(genome.chr[i].getvReadref(Strand::REV), chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);  
  }
};

class shiftJacBit : public ReadShiftProfileGenome {
 public:
 shiftJacBit(const SSPstats &sspst, const SeqStatsGenome &genome):
  ReadShiftProfileGenome("Jaccard index", sspst, genome) {}

  void setDist(ReadShiftProfile &chr, boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev);
  void execchr(const SeqStatsGenome &genome, int32_t i) {
    auto fwd = genBitset(genome.chr[i].getvReadref(Strand::FWD), chr[i].start, chr[i].end);
    auto rev = genBitset(genome.chr[i].getvReadref(Strand::REV), chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftCcp : public ReadShiftProfileGenome {
 public:
 shiftCcp(const SSPstats &sspst, const SeqStatsGenome &genome):
  ReadShiftProfileGenome("Cross correlation", sspst, genome) {}

  void setDist(ReadShiftProfile &chr, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev);
  void execchr(const SeqStatsGenome &genome, int32_t i) {
    auto fwd = genVector(genome.chr[i].getvReadref(Strand::FWD),  chr[i].start, chr[i].end);
    auto rev = genVector(genome.chr[i].getvReadref(Strand::REV),  chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftHamming : public ReadShiftProfileGenome {
 public:
 shiftHamming(const SSPstats &sspst, const SeqStatsGenome &genome):
  ReadShiftProfileGenome("Hamming distance", sspst, genome) {}

  void setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev);
  void execchr(const SeqStatsGenome &genome, int32_t i) {
    auto fwd = genBitset(genome.chr[i].getvReadref(Strand::FWD),  chr[i].start, chr[i].end);
    auto rev = genBitset(genome.chr[i].getvReadref(Strand::REV), chr[i].start, chr[i].end);
    
    setDist(chr[i], fwd, rev);
  }
};

#endif /* _SSP_SHIFTPROFILE_P_H_ */
