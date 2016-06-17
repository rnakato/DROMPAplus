/* 
 * This file was obtained from Q and modified.
 * Copyright (c) 2015, Peter Robinson and Peter Hansen - Charite Universitätsmedizin Berlin
 * http://compbio.charite.de/
 */
#include "pw_hamming.h"
#include "macro.h"
#include "alglib.h"
#include "statistics.h"
#include <map>
#include <boost/dynamic_bitset.hpp>
#include <omp.h>
#include <math.h>   

template <class T>
void GaussianSmoothing(vector<T> &hd)
{
  int size = hd.size();

  vector<double> w(4,0);
  double var=1;
  for(int j=0; j<4; ++j) w[j] = exp((double)(-j*j)/2*var*var);
  double r = 1/(w[0] + (w[1]+w[2]+w[3]) *2);

  double m0;
  double m1(hd[0]);
  double m2(hd[1]);
  double m3(hd[2]);
  for(int i=3; i<size-3; ++i) {
    m0 = hd[i];
    hd[i] = (w[0]*m0 + w[1]*(m1 + hd[i+1]) + w[2]*(m2 + hd[i+2]) + w[3]*(m3 + hd[i+3]))*r;
    m3 = m2;
    m2 = m1;
    m1 = m0;
  }
  
  return;
}

void hammingDistChr(SeqStats &chr, vector<int> &hd, int numthreads)
{
  cout << chr.name <<".." << flush;
  boost::dynamic_bitset<> fwd(chr.len + HD_FROM);
  boost::dynamic_bitset<> rev(chr.len + HD_FROM);
  
  for(int strand=0; strand<STRANDNUM; ++strand) {
#pragma omp parallel for num_threads(numthreads)
    for(uint i=0; i<chr.seq[strand].vRead.size(); ++i) {
      if(chr.seq[strand].vRead[i].duplicate) continue;
      int pos(chr.len -1 - chr.seq[strand].vRead[i].F3);
      if(!RANGE(pos, 0, chr.len-1)) continue;
      if(strand==STRAND_PLUS) fwd.set(pos + HD_FROM);
      else                    rev.set(pos);
    }
  }
    /*for (auto x: chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      int pos(chr.len -1 -x.F3);
      if(!RANGE(pos, 0, chr.len-1)) continue;
      if(strand==STRAND_PLUS) fwd.set(pos + HD_FROM);
      else                    rev.set(pos);
      }
      }*/
  for(int i=0; i<HD_WIDTH; ++i) {
    (fwd >>= 1);
    hd[i]=((fwd ^ rev).count());
  }
  return;
}

void hammingDist(Mapfile &p, int numthreads)
{
  cout << "Estimate fragment length.." << flush;

  if(p.lchr->len > NUM_10M) hammingDistChr(*p.lchr, p.dist.hd, numthreads);
  else {
    for(auto &chr: p.chr) hammingDistChr(chr, p.dist.hd, numthreads);
  }
  
  GaussianSmoothing(p.dist.hd);
  
   // get fragment length FL and HD[FL] run through from (i_num-1),...,2*read_len+1
  int min_hd_fl=p.dist.hd[HD_WIDTH-1];
  int max_hd_fl=p.dist.hd[HD_WIDTH-1];
  int bd = 2*p.dist.lenF3;
  for(int i=HD_WIDTH-1; i>=bd+1; --i) {
    if(p.dist.hd[i] < min_hd_fl) {
      min_hd_fl = p.dist.hd[i];
      p.dist.eflen = i - HD_FROM;
    }
    if(max_hd_fl < p.dist.hd[i]) max_hd_fl = p.dist.hd[i];
  }
  
  // get phantom peak RL and HD[RL] run through from 1,...,2*read_len
  /*  int rl=0;
  int min_hd_rl=p.dist.hd[0];
  for(int i=1; i<=bd; ++i) {
    if(p.dist.hd[i] < min_hd_rl) {
      min_hd_rl = p.dist.hd[i];
      rl = i;
    }
    }*/

  long sum = accumulate(p.dist.hd.begin(), p.dist.hd.end(), 0);
  string filename = p.oprefix + ".hdp.csv";
  ofstream out(filename);
  out << "Strand shift\tHamming distance\tProportion" << endl;
  for(int i=0; i<HD_WIDTH; ++i) out << (i - HD_FROM) << "\t" << p.dist.hd[i] << "\t" << p.dist.hd[i]/(double)sum << endl;
  
  //  double RSC=(max_hd_fl-p.dist.hd[fl])/(double)(max_hd_fl-p.dist.hd[rl]);
  //  out << "RSC: " << RSC << endl;

  cout << "done." << endl;
  cout << "Estimated fragment length: " << p.dist.eflen << endl;
  
  return;
}

void func(short *x, int num, double &ave, double &var)
{
  double dx;
  ave=0; var=0;
  
  for(int i=0; i<num; ++i) ave += x[i];
  ave /= (double)num;
  for(int i=0; i<num; ++i) {
    dx = x[i] - ave;
    var += dx * dx;
  }
  var = sqrt(var);
}

void pw_ccp(Mapfile &p, int numthreads)
{
  printf("Making cross-correlation profile...\n");
  short *fwd = (short *)calloc(p.lchr->len, sizeof(short));
  short *rev = (short *)calloc(p.lchr->len, sizeof(short));

  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto x: p.lchr->seq[strand].vRead) {
      if(x.duplicate) continue;
      int pos(p.lchr->len-1-x.F3);
      if(!RANGE(pos, 0, p.lchr->len-1)) continue;
      if(strand==STRAND_PLUS) ++fwd[pos];
      else                    ++rev[pos];
    }
  }

  int start = HD_FROM;
  int num = p.lchr->len-HD_WIDTH;

  double mx,xx;
  double my,yy;
  func(fwd, num, mx, xx);
  func(rev, num, my, yy);
  map<int, double> mp;
  
#pragma omp parallel for num_threads(numthreads)
  for(int i=-HD_FROM; i<HD_WIDTH; i+=5) {
    double xy(0);
    for(int j=0; j<num; ++j) xy += (fwd[j +start +i] - mx) * (rev[j +start] - my);
    mp[i] = xy / (xx*yy);
  }

  string filename = p.oprefix + ".ccp.csv";
  ofstream out(filename);
  out << "Strand shift\tCross correlation" << endl;
  for(auto itr = mp.begin(); itr != mp.end(); ++itr) out << itr->first << "\t" << itr->second << endl;
  
  free(fwd);
  free(rev);
  return;
}
