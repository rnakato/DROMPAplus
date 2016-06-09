/* 
 * This file was obtained from Q and modified.
 * Copyright (c) 2015, Peter Robinson and Peter Hansen - Charite Universitätsmedizin Berlin
 * http://compbio.charite.de/
 * All rights reserved.
 */

#include "pw_hamming.h"
#include "macro.h"
#include <boost/dynamic_bitset.hpp>

void hammingDist(Mapfile &p)
{
  for(auto chr: p.chr) {
    cout << chr.name << ", "<< chr.len << endl;
    boost::dynamic_bitset<> fwd(chr.len);
    boost::dynamic_bitset<> rev(chr.len);

    for(int strand=0; strand<STRANDNUM; ++strand) {
      for (auto x:chr.seq[strand].vRead) {
	if(x.duplicate) continue;
	int pos(chr.len-1-x.F3);
	if(!RANGE(pos, 0, chr.len-1)) continue;
	if(strand==STRAND_PLUS) fwd.set(pos);
	else                    rev.set(pos);
      }
    }

    map<int, int> hd;
    for(int step=0;step<1000;++step) {
      (fwd >>= 1);
      hd[step]=((fwd ^ rev).count());
    }

    string filename = p.oprefix + ".hammingdistance.xls";
    ofstream out(filename);
    for(int i=0;i<1000;++i) out << i << "\t" << hd[i] << endl;
    break; // zantei
  }

  return;
}
