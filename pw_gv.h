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

class seq{
  long n_read_infile;     /* allow multiread (-kn) */
  long double n_readname; /* number of unique readname */
  long n_read_nonred;     /* number of nonredundant reads (use as "reads number") */
  long n_read_rpkm;       /* number of reads after normalization (n_read_nonred * rpkm weight) */
  long n_read_red;        /* number of redundant reads (filtered as PCR bias) */
  double n_read_afterGC;  /* number of reads after GC normalization */
};

class SeqStats{
  seq single[STRANDNUM];
  seq both;
  double depth;
  double w;
  
  /* FRiP */
  long n_read_inbed;
  double FRiP;
};

class CompStats{
  int nt_all, nt_nonred, nt_red; 
  double complexity;
  int tv;
};

typedef struct{
  SeqStats *genome, *chr;
  Readarray **readarray;
  FragStats fstats;
  WigStats wstats;
  CompStats cs_raw, cs_nonred;
  int threshold4filtering;

  /* for GC*/
  int *GCdist;
  int maxGC;
  int sum_GCdist;
  double *GCweight;
} Mapfile;

#endif /* _PW_GV_H_ */
