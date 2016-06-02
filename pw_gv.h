/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_GV_H_
#define _PW_GV_H_

#include "seq.h"
#include <seqan/bam_io.h>

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

class Readarray{
  int *F3;
  int *F5;
  int *weight;
  //  bool *delete; // for filtering redundant reads
  //bool *ignore; // for ignoring peak regions

  int narray;
};

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

class CompStats{
  int nt_all, nt_nonred, nt_red; 
  double complexity;
  int tv;
};

/*typedef struct{
  SeqStats *genome, *chr;
  Readarray **readarray;
  FragStats fstats;
  WigStats wstats;
  CompStats cs_raw, cs_nonred;
  int threshold4filtering;

  int *GCdist;
  int maxGC;
  int sum_GCdist;
  double *GCweight;
} Mapfile;*/

class Dist{
 public:
  vector<int> readlen;
  vector<int> fraglen; /* +1: over DIST_FRAGLEN_MAX */
  Dist() {
    vector<int> v(DIST_READLEN_MAX,0);
    readlen = v;
    vector<int> v2(DIST_FRAGLEN_MAX,0);
    fraglen = v2;
  }
  void print() {
    int i;
    for(i=0; i<DIST_READLEN_MAX; ++i) if(readlen[i]) cout << i << "\t" << readlen[i]<< endl;
    for(i=0; i<DIST_FRAGLEN_MAX; ++i) if(fraglen[i]) cout << i << "\t" << fraglen[i]<< endl;
  }
};

class FragmentSingle {
public:
  seqan::CharString name;
  int chr;
  int F3;
  Strand strand;
  int fraglen;
  int readlen_F3;

 FragmentSingle(vector<string> v, int flen=0)
  {
    /*    sv = atoi(v[1]); // bitwise FLAG
#ifdef READSV
    if(check_sv(sv, p, clm[0].str)) continue;
#else
    if(check_sv(sv, p)) continue;
#endif
string chr = addchr(v[2]);*/
  }
 FragmentSingle(const seqan::BamAlignmentRecord &record, int flen=0):
  name(record.qName),
    chr(record.rID),
    fraglen(flen),
    readlen_F3(length(record.seq))
  {
    if(!fraglen) fraglen = abs(record.tLen);
    if(hasFlagRC(record)) {
      strand = STRAND_MINUS;
      F3 = record.beginPos + length(record.seq) -1;
    } else {
      strand = STRAND_PLUS;
      F3 = record.beginPos;
    }
  }
  void print() {
    //    cout << name << ", " << chr_F3 << ", " << F3<< ", "<< strand << ", " << "fraglen " << fraglen << "," <<readlen_F3 << endl;
    cout << chr << ", " << F3<< ", "<< strand << ", " << "fraglen " << fraglen << "," <<readlen_F3 << endl;
  }
};

class FragmentPair: public FragmentSingle {
public:
  int chr_F5;
  int F5;
  int readlen_F5;
  //  short num_multimapped;
 FragmentPair(const seqan::BamAlignmentRecord &record):
    FragmentSingle(record)
  {}
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
};

class SeqStats {
 public:
  strandData seq[STRANDNUM];
  double depth;
  double w;
  /* FRiP */
  long nread_inbed;
  double FRiP;
  
 SeqStats(): depth(0), w(0), nread_inbed(0), FRiP(0) {}
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
};
  
class Mapfile {
public:
  Dist dist;
  SeqStats genome;
  vector<SeqStats> chr;

  Mapfile(const RefGenome &g) {
    vector<SeqStats> v(g.chrnum);
    chr = v;
  }
  void addfrag(const FragmentSingle &frag) {
    dist.readlen[frag.readlen_F3]++;
    dist.fraglen[frag.fraglen]++;
    chr[frag.chr].addfrag(frag);
  }
  void update() {
    for (auto x:chr) {
      genome.add(x);
    }
  }
  long nread() {
    long n(0);
    for (auto x:chr) n += x.bothnread();
    return n;
  }
};

#endif /* _PW_GV_H_ */
