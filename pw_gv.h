/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_GV_H_
#define _PW_GV_H_

#include <fstream>
#include <boost/format.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/thread.hpp>
#include "seq.h"
#include "common.h"
#include "util.h"
#include "warn.h"
#include "macro.h"
#include "readdata.h"
#include "statistics.h"

using namespace std;
using namespace boost::program_options;

/* default parameter */
#define NUM_GCOV 5000000

class SeqStats;
void calcFRiP(SeqStats &, const vector<bed>);

class Dist{
  enum {ReadMax=200, 
	FragMax=1000};
 public:
  int lenF3;
  int lenF5;
  int eflen;
  vector<int> readlen;
  vector<int> readlen_F5;
  vector<int> fraglen;

 Dist(): lenF3(0), lenF5(0), eflen(0),
    readlen(ReadMax,0),
    readlen_F5(ReadMax,0),
    fraglen(FragMax,0) {}
  void setlenF3() { lenF3 = getmaxi(readlen); }
  void setlenF5() { lenF5 = getmaxi(readlen_F5); }
  void setFraglen() { eflen = getmaxi(fraglen); }
};

class Fragment {
public:
  string chr;
  int F3;
  Strand strand;
  int fraglen;
  int readlen_F3;

 Fragment(): fraglen(0), readlen_F3(0) {}
 void addSAM(const vector<string> &v, const bool pair, const int sv) {
   chr = rmchr(v[2]);
   readlen_F3 = v[9].length();
   if(pair) fraglen = abs(stoi(v[8]));
   if(sv&16) {
     strand = STRAND_MINUS;
     F3 = stoi(v[3]) + readlen_F3 -1; //SAM
   } else {
     strand = STRAND_PLUS;
     F3 = stoi(v[3]) -1;  //SAM
   }
 }
 void print() const {
   BPRINT("chr:%1%\tposi:%2%\tstrand:%3%\tfraglen:%4%\treadlen:%5%\n") % chr % F3 % strand % fraglen % readlen_F3;
  }
};

class Read {
  int weight;
  enum {WeightNum=1000};
 public:
  int F3;
  int F5;
  int duplicate;
  int inpeak;
  
 Read(const Fragment &frag): weight(WeightNum), F3(frag.F3), duplicate(0), inpeak(0) {
    if(frag.strand == STRAND_PLUS) F5 = frag.F3 + frag.fraglen;
    else F5 = frag.F3 - frag.fraglen;
  }
  double getWeight() const {
    return weight/(double)WeightNum;
  }
  void multiplyWeight(const double w) {
    weight *= w;
  }
};

class strandData {
 public:
  vector<Read> vRead;
  long nread;
  long nread_nonred;
  long nread_red;
  double nread_rpm;
  double nread_afterGC;

 strandData(): nread(0), nread_nonred(0), nread_red(0), nread_rpm(0), nread_afterGC(0) {}
  void setnread() { nread = vRead.size(); }
  void print() const {
    cout << nread << "\t" << nread_nonred << "\t" << nread_red << "\t" << nread_rpm << "\t" << nread_afterGC << endl;
  }
  void printnonred(ofstream &out)  const { printr(out, nread_nonred,  nread); }
  void printred(ofstream &out)     const { printr(out, nread_red,     nread); }
  void printafterGC(ofstream &out) const { printr(out, nread_afterGC, nread); }
  void addReadAfterGC(const double w, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    nread_afterGC += w;
  }
};

// statistics for wigarray
class Wigstats {
  enum{n_mpDist=20, n_wigDist=200};
  long sum;
 public:
  double ave, var, nb_p, nb_n, nb_p0; // pois_p0
  vector<long> mpDist;
  vector<long> wigDist;
  vector<double> pwigDist;

 Wigstats(): sum(0),ave(0), var(0), nb_p(0), nb_n(0), nb_p0(0),
    mpDist(n_mpDist,0), wigDist(n_wigDist,0), pwigDist(n_wigDist,0) {}

  double getPoisson(const int i) const {
    return _getPoisson(i, ave);
  }
  double getNegativeBinomial(const int i) const {
    return _getNegativeBinomial(i, nb_p, nb_n);
  }
  /*  double getZIP(const int i) const {
    return _getZIP(i, ave, pois_p0);
    }*/
  double getZINB(const int i) const {
    return _getZINB(i, nb_p, nb_n, nb_p0);
  }
  int getwigDistthre() const {
    int thre(9);
    long num;
    do{
      ++thre;
      num=0;
      for(int i=0; i<thre; ++i) num += wigDist[i];
    } while(num < sum*0.8);
#ifdef DEBUG
    BPRINT("\nthre %1%  (%2% / %3%)\n") % thre % num % wigDist.size();
#endif
    return thre;
  }
  void estimateParam() {
    int thre = getwigDistthre();
    double par[thre+1];
    par[0] = thre;
    for(int i=0; i<thre; ++i) par[i+1] = pwigDist[i];
    iterateZINB(&par, nb_p, nb_n, nb_p, nb_n, nb_p0);
  }
  void setpWigDist() {
    for(size_t i=0; i<wigDist.size(); ++i) pwigDist[i] = wigDist[i]/(double)sum;
  }
  void getWigStats(const vector<int> &wigarray) {
    int num95 = getPercentile(wigarray, 0.95);
    
    int size = wigDist.size();
    vector<int> ar;
    for(auto x: wigarray) {
      ++sum;
      int v = WIGARRAY2VALUE(x);
      if(v < size) ++wigDist[v];
      if(x >= num95) continue;
      ar.push_back(v);
    }
    setpWigDist();

    moment<int> mm(ar, 0);
    ave = mm.getmean();
    var = mm.getvar();
    nb_p = ave/var;
    if(nb_p>=1) nb_p = 0.9;
    if(nb_p<=0) nb_p = 0.1; 
    nb_n = ave * nb_p /(1 - nb_p);

    estimateParam();
  }
  void printwigDist(ofstream &out, const int i) const {
    out << boost::format("%1%\t%2%\t") % wigDist[i] % pwigDist[i];
  }
  void addmpDist(const double p) {
    if(!RANGE(p,0,1)) cout << "Warning: mappability " << p << " should be [0,1]" << endl;
    else ++mpDist[(int)(p*n_mpDist)];
  }
  void addWigDist(const Wigstats &x) {
    for(uint i=0; i<wigDist.size(); ++i) wigDist[i] += x.wigDist[i];
    sum += x.sum;
    setpWigDist();
  }
  void printmpDist() const {
    long num = accumulate(mpDist.begin(), mpDist.end(), 0);
    for(size_t i=0; i<mpDist.size(); ++i)
      BPRINT("~%1%%%\t%2%\t%3%\n") % ((i+1)*100/mpDist.size()) % mpDist[i] % (mpDist[i]/(double)num); 
  }
  void printPoispar(ofstream &out) const {
    out << boost::format("%1$.3f\t%2$.3f\t") % ave % var;
  }
  void printZINBpar(ofstream &out) const {
    out << boost::format("%1%\t%2%\t%3%") % nb_p % nb_n % nb_p0;
  }
};

class SeqStats {
  bool yeast;
 public:
  string name;
  long len, len_mpbl;
  int nbin;
  double p_mpbl;  /* mappability */
  // genome coverage
  long nbp, ncov, ncovnorm;
  double gcovRaw, gcovNorm;

  Wigstats ws;
  
  strandData seq[STRANDNUM];
  double depth;
  double w;
  /* FRiP */
  long nread_inbed;
  double FRiP;
  
 SeqStats(string s, int l=0): yeast(false), len(l), len_mpbl(l), nbin(0), p_mpbl(0), nbp(0), ncov(0), ncovnorm(0), gcovRaw(0), gcovNorm(0), depth(0), w(0), nread_inbed(0), FRiP(0) {
    name = rmchr(s);
  }
  void addfrag(const Fragment &frag) {
    Read r(frag);
    seq[frag.strand].vRead.push_back(r);
  }
  long bothnread ()         const { return seq[STRAND_PLUS].nread         + seq[STRAND_MINUS].nread; }
  long bothnread_nonred ()  const { return seq[STRAND_PLUS].nread_nonred  + seq[STRAND_MINUS].nread_nonred; }
  long bothnread_red ()     const { return seq[STRAND_PLUS].nread_red     + seq[STRAND_MINUS].nread_red; }
  long bothnread_rpm ()     const { return seq[STRAND_PLUS].nread_rpm     + seq[STRAND_MINUS].nread_rpm; }
  long bothnread_afterGC () const { return seq[STRAND_PLUS].nread_afterGC + seq[STRAND_MINUS].nread_afterGC; }

  void addnread(const SeqStats &x) { 
    for(int i=0; i<STRANDNUM; i++) seq[i].nread += x.seq[i].nread;
  }
  void addnread_red(const SeqStats &x) { 
    for(int i=0; i<STRANDNUM; i++) {
      seq[i].nread_nonred += x.seq[i].nread_nonred;
      seq[i].nread_red    += x.seq[i].nread_red;
    }
  }
  void addGcov(const SeqStats &x, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    nbp      += x.nbp;
    ncov     += x.ncov;
    ncovnorm += x.ncovnorm;
    gcovRaw  = nbp ? ncov / (double)nbp: 0;
    gcovNorm = nbp ? ncovnorm / (double)nbp: 0;
  }

  void print() const {
    cout << name << "\t" << len << "\t" << len_mpbl << "\t" << bothnread() << "\t" << bothnread_nonred() << "\t" << bothnread_red() << "\t" << bothnread_rpm() << "\t" << bothnread_afterGC()<< "\t" << depth << endl;
  }
  void calcdepth(const int flen) {
    depth = len_mpbl ? bothnread_nonred() * flen / (double)len_mpbl: 0;
  }
  void calcGcov(const vector<char> &array) {
    for(long i=0; i<len; ++i) {
      if(array[i] >= MAPPABLE)     ++nbp;      // MAPPABLE || COVREAD_ALL || COVREAD_NORM
      if(array[i] >= COVREAD_ALL)  ++ncov;     // COVREAD_ALL || COVREAD_NORM
      if(array[i] == COVREAD_NORM) ++ncovnorm;
    }
    gcovRaw  = nbp ? ncov     / (double)nbp: 0;
    gcovNorm = nbp ? ncovnorm / (double)nbp: 0;
    //    cout << nbp << ","<< ncov << ","<< ncovnorm << ","<< gcovRaw << ","<< gcovNorm << endl;
  }

  void setF5(int flen) {
    int d;
    for(int strand=0; strand<STRANDNUM; ++strand) {
      if(strand == STRAND_PLUS) d = flen; else d = -flen;
      for(auto &x: seq[strand].vRead) x.F5 = x.F3 + d;
    }
  }
  void setFRiP() {
    FRiP = nread_inbed/(double)bothnread_nonred();
  }
  void setWeight(double weight) {
    w = weight;
    for(int i=0; i<STRANDNUM; i++) seq[i].nread_rpm = seq[i].nread_nonred * w;
  }
  
  void yeaston() { yeast = true; }

  bool isautosome() const {
    int chrnum(0);
    try {
      chrnum = stoi(name);
    } catch (std::invalid_argument e) {  // 数値以外
      if(yeast) { 
	if(name=="I")         chrnum = 1;
	else if(name=="II")   chrnum = 2;
	else if(name=="III")  chrnum = 3;
	else if(name=="IV")   chrnum = 4;
	else if(name=="V")    chrnum = 5;
	else if(name=="VI")   chrnum = 6;
	else if(name=="VII")  chrnum = 7;
	else if(name=="VIII") chrnum = 8;
	else if(name=="IX")   chrnum = 9;
	else if(name=="X")    chrnum = 10;
	else if(name=="XI")   chrnum = 11;
	else if(name=="XII")  chrnum = 12;
	else if(name=="XIII") chrnum = 13;
	else if(name=="XIV")  chrnum = 14;
	else if(name=="XV")   chrnum = 15;
	else if(name=="XVI")  chrnum = 16;
      }
      if(name=="2L") chrnum = 1;
      if(name=="2R") chrnum = 2;
      if(name=="3L") chrnum = 3;
      if(name=="3R") chrnum = 4;
    }
    if(chrnum) return true;
    else       return false;
  }
};

class sepchr {
 public:
  uint s;
  uint e;
 sepchr(uint start, uint end): s(start), e(end) {}
};

class Mapfile {
  bool yeast;
 public:
  string oprefix;
  Dist dist;
  SeqStats genome;
  vector<SeqStats> chr;
  vector<SeqStats>::iterator lchr; // longest chromosome
  vector<sepchr> vsepchr;

  int flen_def;

  // PCR bias
  int thre4filtering;
  int nt_all, nt_nonred, nt_red;
  int tv, gv;
  double r4cmp;
  vector<bed> vbed;
  vector<Peak> vPeak;

  // GC bias
  vector<double> GCweight;
  int maxGC;

  // Wigdist
  int nwigdist;
  vector<int> wigDist;

  Mapfile(const variables_map &values);
  void readGenomeTable(const variables_map &values);
  void getMpbl(const variables_map &values);
  vector<sepchr> getVsepchr(const int);
  
  void addF5(const int readlen_F5) { ++dist.readlen_F5[readlen_F5]; }
  void addfrag(const Fragment &frag) {
    ++dist.readlen[frag.readlen_F3];
    ++dist.fraglen[frag.fraglen];
    int on(0);
    for(auto &x:chr) {
      if(x.name == frag.chr) {
	x.addfrag(frag);
	on++;
      }
    }
    if(!on) cerr << "Warning: " << frag.chr << " is not in genometable." << endl;
  }
  void calcdepth(const variables_map &values) {
    int flen(getflen(values));
    for (auto &x:chr) x.calcdepth(flen);
    genome.calcdepth(flen);
  }
  void setF5(const variables_map &values) {
    int flen(getflen(values));
    for (auto &x:chr) x.setF5(flen);
  }
  void setnread() {
    for (auto &x:chr) {
      for(int i=0; i<STRANDNUM; i++) x.seq[i].setnread();
      genome.addnread(x);
    }
  }
  void setnread_red() {
    for (auto &x:chr) genome.addnread_red(x);
  }
  double complexity() const { return nt_nonred/(double)nt_all; }
  void printstats() const {
    cout << "name\tlength\tlen_mpbl\tread num\tnonred num\tred num\tnormed\tafterGC\tdepth" << endl;
    genome.print();
    for (auto x:chr) x.print();
  }
  void setFRiP() {
    cout << "calculate FRiP score.." << flush;
    for(auto &c: chr) {
      calcFRiP(c, vbed);
      genome.nread_inbed += c.nread_inbed;
    }
    genome.FRiP = genome.nread_inbed/(double)genome.bothnread_nonred();
    
    cout << "done." << endl;
    return;
  }
  int getflen(const variables_map &values) const {
    int flen;
    if(!values.count("nomodel") || values.count("pair")) flen = dist.eflen;
    else flen = flen_def;
    return flen;
  }

  void estimateZINB() {
    int thre = genome.ws.getwigDistthre();
    double par[thre+1];
    par[0] = thre;
    for(int i=0; i<thre; ++i) par[i+1] = genome.ws.pwigDist[i];

    //    iteratePoisson(&par, lchr->ave, genome.ave, genome.pois_p0);
    iterateZINB(&par, lchr->ws.nb_p, lchr->ws.nb_n, genome.ws.nb_p, genome.ws.nb_n, genome.ws.nb_p0);

    return;
  }
};

#endif /* _PW_GV_H_ */
