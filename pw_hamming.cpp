/* 
 * This file was obtained from Q and modified.
 * Copyright (c) 2015, Peter Robinson and Peter Hansen - Charite Universitätsmedizin Berlin
 * http://compbio.charite.de/
 */

#include "pw_hamming.h"
#include "macro.h"
#include "statistics.h"
#include <boost/dynamic_bitset.hpp>

void hammingDist(Mapfile &p)
{
  cout << "Estimate fragment length.." << flush;
  
  //  for(auto &chr: p.chr) {
  //  cout << p.lchr->name <<".." << flush;
  boost::dynamic_bitset<> fwd(p.lchr->len + HD_FROM);
  boost::dynamic_bitset<> rev(p.lchr->len + HD_FROM);

  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto x: p.lchr->seq[strand].vRead) {
      if(x.duplicate) continue;
      int pos(p.lchr->len-1-x.F3);
      if(!RANGE(pos, 0, p.lchr->len-1)) continue;
      if(strand==STRAND_PLUS) fwd.set(pos + HD_FROM);
      else                    rev.set(pos);
    }
  }
  
  for(int i=0; i<HD_WIDTH; ++i) {
    (fwd >>= 1);
    p.dist.hd[i]=((fwd ^ rev).count());
  }
  //  }
  
  //  vector<int> hdg(HD_WIDTH,0);
  //for(int i=0; i<HD_WIDTH; ++i) {
    //    for(auto chr: p.chr) hdg[i] += chr.hd[i];
  //   hdg[i] = p.dist.hd[i];
  // }

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

void pw_ccp(Mapfile &p)
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

  /*
  TYPE_WIGARRAY qnt99 = calc_qnt(plus, chrlen, 0.99);
  for(i=0; i<chrlen; i++){
    if(plus[i]  > qnt99) plus[i]  = qnt99;
    if(minus[i] > qnt99) minus[i] = qnt99;
    }*/

  string filename = p.oprefix + ".ccp.csv";
  ofstream out(filename);
  out << "Strand shift\tCross correlation\tProportion" << endl;
  int start = HD_FROM;
  int num = p.lchr->len-HD_WIDTH;

  double mx,xx;
  double my,yy;
  func(fwd, num, mx, xx);
  func(rev, num, my, yy);

  for(int i=-HD_FROM; i<HD_WIDTH; i+=5) {
    double xy(0);
    for(int j=0; j<num; ++j) xy += (fwd[j +start +i] - mx) * (rev[j +start] - my);
    out << i << "\t" << xy / (xx*yy) << endl;
  }

  free(fwd);
  free(rev);
  return;
}
