/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_GV_H_
#define _PW_GV_H_

#include "seq.h"
#include "common.h"
#include "util.h"
//#include <seqan/bam_io.h>
#include <map>
#include <fstream>
#include <boost/format.hpp>

/* default parameter */
#define FRAGMENT_LEN 150
#define MAX_FRAGMENT_LEN 500
#define NUM4RPM_DEFAULT 20000000
#define NUM4DEPTH_DEFAULT 1.0
#define NUM4CMP_DEFAULT 10
#define FLEN4GC_DEFAULT 120

#define DIST_READLEN_MAX 200
#define DIST_FRAGLEN_MAX 1000

#define NUM_DARRAY 100
#define READARRAY_NUM 50000 

class WigStatsMember {
  int max;
  int *darray_all, *darray_bg;
  int num;
  double ave, var;
  double nb_p, nb_n, nb_p0;  /* for negative binomial model */
};

class WigStats{
  int n_darray;
  int thre;
  int num95;
  WigStatsMember *genome, *chr;
};

class Dist{
 public:
  vector<int> readlen;
  vector<int> readlen_F5;
  vector<int> fraglen; /* +1: over DIST_FRAGLEN_MAX */
  int eflen;
 Dist(): eflen(0) {
    vector<int> v(DIST_READLEN_MAX,0);
    readlen = v;
    vector<int> v2(DIST_READLEN_MAX,0);
    readlen_F5 = v2;
    vector<int> v3(DIST_FRAGLEN_MAX,0);
    fraglen = v3;
  }
  void printReadlen(string &outputfile, bool pair) {
    long sum(0);
    for(int i=0; i<DIST_READLEN_MAX; ++i) sum += readlen[i];
    ofstream out(outputfile);
    out << "F3 read length distribution" << endl;
    out << "length\tread number\tproportion" << endl;
    for(int i=0; i<DIST_READLEN_MAX; ++i) if(readlen[i]) out << boost::format("%1%\t%2%\t%3%\n") % i % readlen[i] % (readlen[i]/(double)sum);

    if(pair) {
      out << "\n\nF5 read length distribution" << endl;
      out << "length\tread number\tproportion" << endl;
      for(int i=0; i<DIST_READLEN_MAX; ++i) if(readlen_F5[i]) out << boost::format("%1%\t%2%\t%3%\n") % i % readlen_F5[i] % (readlen_F5[i]/(double)sum);
    }
  }

  void printFraglen(string &outputfile) {
    int flen_max=0;
    long sum(0);
    for(int i=0; i<DIST_FRAGLEN_MAX; ++i) {
      if(flen_max < fraglen[i]){
	flen_max = fraglen[i];
	eflen = i;
      }
      sum += readlen[i];
    }
    ofstream out(outputfile);
    out << "length\tread number\tproportion" << endl;
    for(int i=0; i<DIST_FRAGLEN_MAX; ++i) if(fraglen[i]) out << boost::format("%1%\t%2%\t%3%\n") % i % fraglen[i] % (fraglen[i]/(double)sum);
  }
};

class FragmentSingle {
public:
  string name;
  string chr;
  int F3;
  Strand strand;
  int fraglen;
  int readlen_F3;

 FragmentSingle(vector<string> v, int flen=0):
  name(v[0]), chr(v[2]), fraglen(flen), readlen_F3(v[9].length())
    {
      if(!fraglen) fraglen = abs(stoi(v[8]));
      int sv = stoi(v[1]); // bitwise FLAG
      if(sv&16) {
	strand = STRAND_MINUS;
	F3 = stoi(v[3]) + readlen_F3 -1; //SAM
      } else {
	strand = STRAND_PLUS;
	F3 = stoi(v[3]) -1;  //SAM
      }
    }
  void print() {
    cout << name << ", " << chr << ", " << F3 << ", "<< strand << ", " << "fraglen " << fraglen << "," <<readlen_F3 << endl;
  }
};

class Read {
 public:
  int F3;
  int F5;
  int weight;
  int duplicate;
  int inpeak;
 Read(const FragmentSingle &frag): F3(frag.F3), weight(1), duplicate(0), inpeak(0) {
    if(frag.strand == STRAND_PLUS) F5 = frag.F3 + frag.fraglen;
    else F5 = frag.F3 - frag.fraglen;
  }
};

class strandData {
 public:
  vector<Read> vRead;
  long nread_nonred;
  long nread_red;
  double nread_rpm;
  double nread_afterGC;
 strandData(): nread_nonred(0), nread_red(0), nread_rpm(0), nread_afterGC(0) {}
  long nread() { return (long)vRead.size();}
  void print() {
    cout << nread() << "\t" << nread_nonred << "\t" << nread_red << "\t" << nread_rpm << "\t" << nread_afterGC << endl;
  }
  void printnonred(ofstream &out) {
    printr(out, nread_nonred, nread());
  }
  void printred(ofstream &out) {
    printr(out, nread_red, nread());
  }
  void printafterGC(ofstream &out) {
    printr(out, nread_afterGC, nread());
  }
};

class SeqStats {
 public:
  string name;
  long len, len_mpbl;
  int nbin;
  double p_mpbl;  /* mappability */
  double gcov;    /* genome coverage for bin */

  strandData seq[STRANDNUM];
  double depth;
  double w;
  /* FRiP */
  long nread_inbed;
  double FRiP;

 SeqStats(string s, int l=0): name(s),len(l), len_mpbl(l), nbin(0), p_mpbl(0), gcov(0), depth(0), w(0), nread_inbed(0), FRiP(0) { }
  void addfrag(const FragmentSingle &frag) {
    Read r(frag);
    seq[frag.strand].vRead.push_back(r);
  }
  long bothnread ()         { return seq[STRAND_PLUS].nread()       + seq[STRAND_MINUS].nread(); }
  long bothnread_nonred ()  { return seq[STRAND_PLUS].nread_nonred  + seq[STRAND_MINUS].nread_nonred; }
  long bothnread_red ()     { return seq[STRAND_PLUS].nread_red     + seq[STRAND_MINUS].nread_red; }
  long bothnread_rpm ()     { return seq[STRAND_PLUS].nread_rpm     + seq[STRAND_MINUS].nread_rpm; }
  long bothnread_afterGC () { return seq[STRAND_PLUS].nread_afterGC + seq[STRAND_MINUS].nread_afterGC; }

  void add(const SeqStats &x) { 
    for(int i=0; i<STRANDNUM; i++) {
      seq[i].nread_nonred  += x.seq[i].nread_nonred;
      seq[i].nread_red     += x.seq[i].nread_red;
      seq[i].nread_rpm     += x.seq[i].nread_rpm;
      seq[i].nread_afterGC += x.seq[i].nread_afterGC;
    }
  }
  void print() {
    cout << name << "\t" << len << "\t" << len_mpbl << "\t" << bothnread() << "\t" << bothnread_nonred() << "\t" << bothnread_red() << "\t" << bothnread_rpm() << "\t" << bothnread_afterGC()<< "\t" << depth << endl;
  }
  void calcdepth(int flen) {
    depth = len_mpbl ? bothnread_nonred() * flen / (double)len_mpbl: 0;
  }
};
  
class Mapfile {
public:
  Dist dist;
  SeqStats genome;
  vector<SeqStats> chr;

  // PCR bias
  int thre4filtering;
  int nt_all, nt_nonred, nt_red;
  int tv;
  double r4cmp;

  Mapfile(string gtfile, int binsize, int flen);
  void addF5(const int readlen_F5) { dist.readlen_F5[readlen_F5]++; }
  void addfrag(const FragmentSingle &frag) {
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
  void calcdepth() {
    for (auto &x:chr) x.calcdepth(dist.eflen);
    genome.calcdepth(dist.eflen);
  }
  void update() { for (auto x:chr) genome.add(x); }
  long nread() {
    long n(0);
    for (auto x:chr) n += x.bothnread();
    return n;
  }
  double complexity() { return nt_nonred/(double)nt_all; }
  void printstats() {
    for (auto x:chr) x.print();
    genome.print();
  }
};

#endif /* _PW_GV_H_ */
