/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_GV_H_
#define _PW_GV_H_

#include "seq.h"

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

class FragStats{
  int dist_readlen_F3[DIST_READLEN_MAX];
  int dist_readlen_F5[DIST_READLEN_MAX];
  int dist_fraglen[DIST_FRAGLEN_MAX +1]; /* +1: over DIST_FRAGLEN_MAX */
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

class strandStats {
public:
  long nread;
  long nread_nonred;
  long nread_red;
  double nread_rpm;
  double nread_afterGC;
  strandStats(): nread(0), nread_nonred(0), nread_red(0), nread_rpm(0), nread_afterGC(0) {}
};

class ChrStats {
public:
  strandStats seq[STRANDNUM];
  strandStats both;
  double depth;
  double w;

  /* FRiP */
  long nread_inbed;
  double FRiP;
 ChrStats(): depth(0), w(0), nread_inbed(0), FRiP(0) {}
};

class SeqStats {
  ChrStats genome;
  vector<ChrStats> chr;
};



#endif /* _PW_GV_H_ */
