/* 
 * This file was obtained from Q and modified.
 * Copyright (c) 2015, Peter Robinson and Peter Hansen - Charite Universitätsmedizin Berlin
 * http://compbio.charite.de/
 */

#include "pw_hamming.h"
#include "macro.h"
#include <boost/dynamic_bitset.hpp>

void hammingDist(Mapfile &p)
{
  cout << "Plot hamming distance.." << flush;
  
  //  for(auto &chr: p.chr) {
  cout << p.lchr->name <<".." << flush;
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
    p.lchr->hd[i]=((fwd ^ rev).count());
  }
  //  }
  
  vector<int> hdg(HD_WIDTH,0);
  for(int i=0; i<HD_WIDTH; ++i) {
    //    for(auto chr: p.chr) hdg[i] += chr.hd[i];
    hdg[i] = p.lchr->hd[i];
  }

  // get fragment length FL and HD[FL] run through from (i_num-1),...,2*read_len+1
  int fl=0;
  int min_hd_fl=hdg[HD_WIDTH-1];
  int max_hd_fl=hdg[HD_WIDTH-1];
  int bd = 2*p.dist.lenF3;
  for(int i=HD_WIDTH-1; i>=bd+1; i--) {
    if(hdg[i]<min_hd_fl) {
      min_hd_fl=hdg[i];
      fl=i;
    }
    if(max_hd_fl<hdg[i]) max_hd_fl=hdg[i];
  }
  
  // get phantom peak RL and HD[RL] run through from 1,...,2*read_len
  int rl=0;
  int min_hd_rl=hdg[0];
  for(int i=1; i<=bd; ++i) {
    if(hdg[i]<min_hd_rl) {
      min_hd_rl=hdg[i];
      rl=i;			
    }
  }

  // get RSC
  double RSC=(max_hd_fl-hdg[fl])/(double)(max_hd_fl-hdg[rl]);

  string filename = p.oprefix + ".hdp.csv";
  ofstream out(filename);
  for(int i=0; i<HD_WIDTH; ++i) out << (i - HD_FROM) << "\t" << hdg[i] << endl;
  out << "RSC: " << RSC << endl;

  cout << "done." << endl;
  
  return;
}
