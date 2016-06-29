/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_GV_H_
#define _PW_GV_H_

#include <fstream>
#include <boost/format.hpp>
#include <omp.h>
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
#define FRAGMENT_LEN 150
#define MAX_FRAGMENT_LEN 500
#define NUM4RPM_DEFAULT 20000000
#define NUM4DEPTH_DEFAULT 1.0
#define NUM4CMP_DEFAULT 10
#define NUM_GCOV 5000000
#define READWEIGHTNUM 1000

#define FLEN4GC_DEFAULT 120
#define DIST_READLEN_MAX 200
#define DIST_FRAGLEN_MAX 1000

#define NUM_WIGDISTARRAY 200
#define NUM_MPDIST 20

#define HD_FROM 200
//#define HD_WIDTH 5000
#define HD_WIDTH 600

class shiftDist{
  double bk;
  
 public:
  map<int, double> mp;
  map<int, double> nc;
  double r;
  double nsc;
  int nsci;
  
 shiftDist(): bk(0), r(0), nsc(0), nsci(0) {}

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
  void getflen(int lenF3) {
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
};

class Dist{
 public:
  int lenF3;
  int lenF5;
  int eflen;
  vector<int> readlen;
  vector<int> readlen_F5;
  vector<int> fraglen;
  vector<int> hd;

 Dist(): lenF3(0), lenF5(0), eflen(0) {
    vector<int> v(DIST_READLEN_MAX,0);
    readlen = v;
    vector<int> v2(DIST_READLEN_MAX,0);
    readlen_F5 = v2;
    vector<int> v3(DIST_FRAGLEN_MAX,0);
    fraglen = v3;
    vector<int> h(HD_WIDTH,0);
    hd = h;
  }
  void setlenF3() { lenF3 = getmaxi(readlen); }
  void setlenF5() { lenF5 = getmaxi(readlen_F5); }
  void setFraglen() { eflen = getmaxi(fraglen); }
};

class Fragment {
public:
  //  string name;
  string chr;
  int F3;
  Strand strand;
  int fraglen;
  int readlen_F3;

 Fragment(): fraglen(0), readlen_F3(0) {}
 void addSAM(const vector<string> &v, const bool pair, const int sv) {
   chr = v[2];
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
 void print() {
   BPRINT("chr:%1%\tposi:%2%\tstrand:%3%\tfraglen:%4%\treadlen:%5%\n") % chr % F3 % strand % fraglen % readlen_F3;
  }
};

class Read {
  int weight;
 public:
  int F3;
  int F5;
  int duplicate;
  int inpeak;
 Read(const Fragment &frag): weight(READWEIGHTNUM), F3(frag.F3), duplicate(0), inpeak(0) {
    if(frag.strand == STRAND_PLUS) F5 = frag.F3 + frag.fraglen;
    else F5 = frag.F3 - frag.fraglen;
  }
  double getWeight() const {
    return weight/(double)READWEIGHTNUM;
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
  void print() {
    cout << nread << "\t" << nread_nonred << "\t" << nread_red << "\t" << nread_rpm << "\t" << nread_afterGC << endl;
  }
  void printnonred(ofstream &out)  const { printr(out, nread_nonred,  nread); }
  void printred(ofstream &out)     const { printr(out, nread_red,     nread); }
  void printafterGC(ofstream &out) const { printr(out, nread_afterGC, nread); }
};

class SeqStats {
  int num95;
 public:
  string name;
  long len, len_mpbl;
  int nbin, nbindist;
  double p_mpbl;  /* mappability */
  // genome coverage
  long nbp, ncov, ncovnorm;
  double gcovRaw, gcovNorm;
  // statistics of wigarray
  vector<long> mpDist;
  vector<long> wigDist;
  double ave, var, nb_p, nb_n, nb_p0, pois_p0;

  strandData seq[STRANDNUM];
  double depth;
  double w;
  /* FRiP */
  long nread_inbed;
  double FRiP;

  shiftDist jac;
  shiftDist ccp;
  double rchr;
  
 SeqStats(string s, int l=0): num95(0), name(s),len(l), len_mpbl(l), nbin(0), nbindist(0), p_mpbl(0), nbp(0), ncov(0), ncovnorm(0), gcovRaw(0), gcovNorm(0), depth(0), w(0), nread_inbed(0), FRiP(0), rchr(0) {
    vector<long> v(NUM_MPDIST,0); // 5% div
    mpDist = v;
    vector<long> v2(NUM_WIGDISTARRAY,0);
    wigDist = v2;
  }
  void addfrag(const Fragment &frag) {
    Read r(frag);
    seq[frag.strand].vRead.push_back(r);
  }
  long bothnread () const         { return seq[STRAND_PLUS].nread         + seq[STRAND_MINUS].nread; }
  long bothnread_nonred () const  { return seq[STRAND_PLUS].nread_nonred  + seq[STRAND_MINUS].nread_nonred; }
  long bothnread_red () const     { return seq[STRAND_PLUS].nread_red     + seq[STRAND_MINUS].nread_red; }
  long bothnread_rpm () const     { return seq[STRAND_PLUS].nread_rpm     + seq[STRAND_MINUS].nread_rpm; }
  long bothnread_afterGC () const { return seq[STRAND_PLUS].nread_afterGC + seq[STRAND_MINUS].nread_afterGC; }

  void addnread(const SeqStats &x) { 
    for(int i=0; i<STRANDNUM; i++) {
      seq[i].nread += x.seq[i].nread;
    }
  }
  void addnread_red(const SeqStats &x) { 
    for(int i=0; i<STRANDNUM; i++) {
      seq[i].nread_nonred += x.seq[i].nread_nonred;
      seq[i].nread_red    += x.seq[i].nread_red;
    }
  }
  void addGcov(const SeqStats &x) {
    nbp      += x.nbp;
    ncov     += x.ncov;
    ncovnorm += x.ncovnorm;
    gcovRaw  = nbp ? ncov / (double)nbp: 0;
    gcovNorm = nbp ? ncovnorm / (double)nbp: 0;
  }

  void addjac(const SeqStats &x) {
    addmp(jac.mp, x.jac.mp, x.rchr);
    addmp(jac.nc, x.jac.nc, x.rchr);
  }
  void addccp(const SeqStats &x) {
    addmp(ccp.mp, x.ccp.mp, x.rchr);
    addmp(ccp.nc, x.ccp.nc, x.rchr);
  }
  void print() {
    cout << name << "\t" << len << "\t" << len_mpbl << "\t" << bothnread() << "\t" << bothnread_nonred() << "\t" << bothnread_red() << "\t" << bothnread_rpm() << "\t" << bothnread_afterGC()<< "\t" << depth << endl;
  }
  void calcdepth(const int flen) {
    depth = len_mpbl ? bothnread_nonred() * flen / (double)len_mpbl: 0;
  }
  void calcGcov(const vector<char> &array) {
    for(int i=0; i<len; ++i) {
      if(array[i] == MAPPABLE    || array[i] == COVREAD_ALL || array[i] == COVREAD_NORM) ++nbp;
      if(array[i] == COVREAD_ALL || array[i] == COVREAD_NORM) ++ncov;
      if(array[i] == COVREAD_NORM) ++ncovnorm;
    }
    gcovRaw  = nbp ? ncov / (double)nbp: 0;
    gcovNorm = nbp ? ncovnorm / (double)nbp: 0;
  }
  void printwigDist(ofstream &out, const int i) {
    out << boost::format("%1%\t%2%\t") % wigDist[i] % (wigDist[i]/(double)nbindist);
  }
  void addmpDist(const double p) {
    if(!RANGE(p,0,1)) cout << "Warning: mappability " << p << " should be [0,1]" << endl;
    else ++mpDist[(int)(p*NUM_MPDIST)];
  }
  void printmpDist() {
    long num =  accumulate(mpDist.begin(), mpDist.end(), 0);
    for(int i=0; i<NUM_MPDIST; ++i) BPRINT("~%1%%%\t%2%\t%3%\n") % ((i+1)*100/NUM_MPDIST) % mpDist[i] % (mpDist[i]/(double)num); 
  }

  void getWigStats(const vector<int> &wigarray) {
    num95 = getPercentile(wigarray, 0.95);
    nbindist = wigarray.size();

    int size=wigDist.size();
    vector<int> ar;
    for(auto x: wigarray) {
      int v = WIGARRAY2VALUE(x);
      if(v < size) ++wigDist[v];
      if(x >= num95) continue;
      ar.push_back(v);
    }
    getMoment(ar, ave, var);
    nb_p = ave/var;
    if(nb_p>=1) nb_p = 0.9;
    if(nb_p<=0) nb_p = 0.1; 
    nb_n = ave * nb_p /(1 - nb_p);
  }
  void addWigStats(const SeqStats &x) {
    nbindist += x.nbindist;
    for(uint i=0; i<wigDist.size(); ++i) wigDist[i] += x.wigDist[i]; 
  }
  void printWigStats() {
    BPRINT("ave=%1%, var=%2%, p=%3%, n=%4%\n") % ave % var % nb_p % nb_n;
  }
  void setF5(int flen) {
    int d;
    for(int strand=0; strand<STRANDNUM; ++strand) {
      if(strand == STRAND_PLUS) d = flen; else d = -flen;
      for (auto &x: seq[strand].vRead) x.F5 = x.F3 + d;
    }
  }
  void setWeight(double weight) {
    w = weight;
    for(int i=0; i<STRANDNUM; i++) seq[i].nread_rpm = seq[i].nread_nonred * w;
  }
  void calcFRiP(const vector<bed> vbed) {
    vector<char> array(len,MAPPABLE);
    arraySetBed(array, name, vbed);
    for(int strand=0; strand<STRANDNUM; ++strand) {
      for (auto &x: seq[strand].vRead) {
	if(x.duplicate) continue;
	int s(min(x.F3, x.F5));
	int e(max(x.F3, x.F5));
	for(int i=s; i<=e; ++i) {
	  if(array[i]==INBED) {
	    x.inpeak = 1;
	    nread_inbed++;
	    break;
	  }
	}
      }
    }
    FRiP = nread_inbed/(double)bothnread_nonred();
  }
  double getPoisson(const int i) {
    return _getPoisson(i, ave);
  }
  double getNegativeBinomial(const int i) {
    return _getNegativeBinomial(i, nb_p, nb_n);
  }
  double getZIP(const int i) {
    return _getZIP(i, ave, pois_p0);
  }
  double getZINB(const int i) {
    return _getZINB(i, nb_p, nb_n, nb_p0);
  }
  int getwigDistthre() {
    int thre(9);
    long num;
    do{
      ++thre;
      num=0;
      for(int i=0; i<thre; ++i) num += wigDist[i];
    } while(num < nbindist*0.8);
#ifdef DEBUG
    BPRINT("\nthre %1%  (%2% / %3%)\n") % thre % num % nbindist;
#endif
    return thre;
  }
  
  void estimateParam() {
    int thre = getwigDistthre();
    double par[thre+1];
    par[0] = thre;
    for(int i=0; i<thre; ++i) par[i+1] = wigDist[i] /(double)nbindist;
    iterateZINB(&par, nb_p, nb_n, nb_p, nb_n, nb_p0);
  }
};

class Mapfile {
public:
  string oprefix;
  Dist dist;
  SeqStats genome;
  vector<SeqStats> chr;
  vector<SeqStats>::iterator lchr; // longest chromosome

  string lastchr;
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

  // WigStats
  int nwigdist;
  vector<int> wigDist;

  Mapfile(const variables_map &values);
  void addF5(const int readlen_F5) { dist.readlen_F5[readlen_F5]++; }
  void addfrag(const Fragment &frag) {
    dist.readlen[frag.readlen_F3]++;
    dist.fraglen[frag.fraglen]++;
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
    for (auto &x:chr) {
      genome.addnread_red(x);
      x.rchr = x.bothnread_nonred()/(double)genome.bothnread_nonred();
    }
  }
  double complexity() const { return nt_nonred/(double)nt_all; }
  void printstats() {
    for (auto x:chr) x.print();
    genome.print();
  }
  void calcFRiP() {
    cout << "calculate FRiP score.." << flush;
    for(auto &c: chr) {
      c.calcFRiP(vbed);
      genome.nread_inbed += c.nread_inbed;
    }
    genome.FRiP = genome.nread_inbed/(double)genome.bothnread_nonred();
    
    cout << "done." << endl;
    return;
  }
  int getflen(const variables_map &values) {
    int flen;
    if(!values.count("nomodel") || values.count("pair")) flen = dist.eflen;
    else flen = flen_def;
    return flen;
  }

  void estimateZINB() {
    int thre = genome.getwigDistthre();
    double par[thre+1];
    par[0] = thre;
    for(int i=0; i<thre; ++i) par[i+1] = genome.wigDist[i] /(double)genome.nbindist;

    //    iteratePoisson(&par, lchr->ave, genome.ave, genome.pois_p0);
    iterateZINB(&par, lchr->nb_p, lchr->nb_n, genome.nb_p, genome.nb_n, genome.nb_p0);

    //    for (auto &x:chr) x.estimateParam(thre);
    
    return;
  }
};

#endif /* _PW_GV_H_ */
