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

#define CHR1ONLY 1
#define HD_NG_FROM 4000
#define HD_NG_TO   5000
#define HD_NG_STEP 100

void jaccardfunc_vector(Mapfile &p, SeqStats &chr);
void jaccardfunc_bitset(Mapfile &p, SeqStats &chr);
void ccpfunc(Mapfile &p, SeqStats &chr);
  
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

void outputmp(ofstream &out, shiftDist &dist, const Mapfile &p, const string str) {
  dist.setControlRatio();
  dist.getflen(p.dist.lenF3);
  double sum(dist.getmpsum());
  
  out << "NSC\t"<< dist.nsc << endl;
  out << "Estimated fragment length\t" << dist.nsci << endl;
  out << "Background enrichment\t" << dist.getBackEnrich(p.genome.bothnread_nonred()) << endl;
  out << "Strand shift\t" << str << "\tprop\tper 10M reads\tper control" << endl;
  for(auto itr = dist.mp.begin(); itr != dist.mp.end(); ++itr) {
    out << itr->first << "\t" << itr->second << "\t" << (itr->second/sum)<< "\t" << (itr->second*NUM_10M/p.genome.bothnread_nonred()) << "\t" << (itr->second * dist.r) << endl;
  }
}

void pw_Jaccard(Mapfile &p, int numthreads)
{
  printf("Making Jaccard index profile...\n");

#ifdef CHR1ONLY
  for(uint i=0; i<1; ++i) {
#else 
#pragma omp parallel for num_threads(numthreads)
  for(uint i=0; i<p.chr.size(); ++i) {
#endif
    jaccardfunc_vector(p, p.chr[i]);
    //jaccardfunc_bitset(p, p.chr[i]);
    
    p.genome.addjac(p.chr[i]);
    string filename = p.oprefix + ".jaccard." + p.chr[i].name +".csv";
    ofstream out(filename);
    outputmp(out, p.chr[i].jac, p, "Jaccard index");
  }

  string filename = p.oprefix + ".jaccard.csv";
  ofstream out(filename);
  outputmp(out, p.genome.jac, p, "Jaccard index");
  
  return;
}
  
void pw_ccp(Mapfile &p, int numthreads)
{
  printf("Making cross-correlation profile...\n");

#ifdef CHR1ONLY
  for(uint i=0; i<1; ++i) {
#else 
#pragma omp parallel for num_threads(numthreads)
  for(uint i=0; i<p.chr.size(); ++i) {
#endif
    ccpfunc(p, p.chr[i]);
  }
  
  string filename = p.oprefix + ".ccp.csv";
  ofstream out(filename);
  outputmp(out, p.genome.ccp, p, "Cross correlation");
  
  return;
}


void jaccardfunc_bitset(Mapfile &p, SeqStats &chr)
{
  cout << chr.name << endl;
  int start(0);
  int end(chr.len);
  //    int start(1.213*NUM_100M);
  //int end(1.214*NUM_100M);
  
  int width(end-start);
  
  boost::dynamic_bitset<> fwd(width + HD_FROM);
  boost::dynamic_bitset<> rev(width + HD_FROM);

  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto x: chr.seq[strand].vRead) {
      if(x.duplicate) continue;

      int pos(chr.len -1 -x.F3);
      if(!RANGE(pos, start, end-1)) continue;
      if(strand==STRAND_PLUS) fwd.set(pos - start + HD_FROM);
      else                    rev.set(pos - start);
    }
  }
  
  int xx = fwd.count();
  int yy = rev.count();
  int xysum(xx+yy);
  
  //#pragma omp parallel for num_threads(numthreads)
  for(int i=-HD_FROM; i<HD_WIDTH; ++i) {
    (fwd >>= 1);
    int xy((fwd & rev).count());
    chr.jac.mp[i] = xy/(double)(xysum-xy);
  }
  
  for(int i=HD_NG_FROM; i<HD_NG_TO; ++i) {
    (fwd >>= 1);
    int xy((fwd & rev).count());
    chr.jac.nc[i] = xy/(double)(xysum-xy);
  }

  return;
}
  
void jaccardfunc_vector(Mapfile &p, SeqStats &chr)
{
  cout << chr.name << endl;
  int start(0);
  int end(chr.len);
  //    int start(1.213*NUM_100M);
  //int end(1.214*NUM_100M);
  
  int width(end-start);
  
  vector<char> fwd(width,0);
  vector<char> rev(width,0);
  
  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto x: chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      if(!RANGE(x.F3, start, end-1)) continue;
      if(strand==STRAND_PLUS) ++fwd[x.F3 - start];
      else                    ++rev[x.F3 - start];
    }
  }
  
  /*    for(int j=HD_FROM; j< width - HD_WIDTH; ++j) {
	if(fwd[j] && rev[j + p.dist.lenF3]) cout << (j+start) << "\t" << (int)fwd[j]<< "\t" <<  (int)rev[j + p.dist.lenF3]<< "\t" << p.dist.lenF3 << endl;
	}*/

  int xx = accumulate(fwd.begin(), fwd.end(), 0);
  int yy = accumulate(rev.begin(), rev.end(), 0);
  int xysum(xx+yy);
  
  //#pragma omp parallel for num_threads(numthreads)
  for(int step=-HD_FROM; step<HD_WIDTH; ++step) {
    int xy(0);
    //#pragma omp parallel for num_threads(numthreads) reduction(+:xy)
    for(int j=HD_FROM; j<width - HD_NG_TO; ++j) if(fwd[j] * rev[j+step]) xy += max(fwd[j], rev[j+step]);
    chr.jac.mp[step] = xy/(double)(xysum-xy);
  }
  
  for(int step=HD_NG_FROM; step<HD_NG_TO; step+=HD_NG_STEP) {
    int xy(0);
    //#pragma omp parallel for num_threads(numthreads) reduction(+:xy)
    for(int j=HD_FROM; j<width - HD_NG_TO; ++j) if(fwd[j] * rev[j+step]) xy += max(fwd[j], rev[j+step]);
    chr.jac.nc[step] = xy/(double)(xysum-xy);
  }

  return;
}
  
template <class T>
void calcMeanSD(const vector<T> &x, int max, double &ave, double &sd)
{
  double dx, var(0);

  ave=0;
  for(int i=HD_FROM; i<max; ++i) ave += x[i];
  ave /= (double)(max - HD_FROM);
  for(int i=HD_FROM; i<max; ++i) {
    dx = x[i] - ave;
    var += dx * dx;
  }
  sd = sqrt(var/double(max -HD_FROM -1));
}


void ccpfunc(Mapfile &p, SeqStats &chr)
{
  cout << chr.name << endl;
  int start(0);
  int end(chr.len);
  //int start(1.5*NUM_100M);
  //int end(2*NUM_100M);
  //    int start(1.213*NUM_100M);
  //int end(1.214*NUM_100M);
  
  int width(end-start);
  
  vector<char> fwd(width,0);
  vector<char> rev(width,0);
  
  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto x: chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      if(!RANGE(x.F3, start, end-1)) continue;
      if(strand==STRAND_PLUS) ++fwd[x.F3 - start];
      else                    ++rev[x.F3 - start];
    }
  }
  
  int max = width - HD_NG_TO;
  double mx,xx;
  double my,yy;
  calcMeanSD(fwd, max, mx, xx);
  calcMeanSD(rev, max, my, yy);

  //#pragma omp parallel for num_threads(numthreads)
  for(int step=-HD_FROM; step<HD_WIDTH; step+=5) {
    double xy(0);
    //#pragma omp parallel for num_threads(numthreads) reduction(+:xy)
    for(int j=HD_FROM; j<max; ++j) xy += (fwd[j] - mx) * (rev[j+step] - my);
    chr.ccp.mp[step] += xy;
  }
  
  for(int step=HD_NG_FROM; step<HD_NG_TO; step+=HD_NG_STEP) {
    double xy(0);
    //#pragma omp parallel for num_threads(numthreads) reduction(+:xy)
    for(int j=HD_FROM; j<max; ++j) xy += (fwd[j] - mx) * (rev[j+step] - my);
    chr.ccp.nc[step] += xy;
  }

  double val(1/(xx*yy*(max - HD_FROM - 1)));
  for(auto itr = chr.ccp.mp.begin(); itr != chr.ccp.mp.end(); ++itr) itr->second *= val;
  for(auto itr = chr.ccp.nc.begin(); itr != chr.ccp.nc.end(); ++itr) itr->second *= val;
  
  p.genome.addccp(chr);

  string filename = p.oprefix + ".ccp." + chr.name +".csv";
  ofstream out(filename);
  outputmp(out, chr.ccp, p, "Cross correlation");

  return;
}
 
void hammingDistChr(SeqStats &chr, vector<int> &hd)
{
  int start(0);
  int end(chr.len);
  int width(end-start);  

  boost::dynamic_bitset<> fwd(width + HD_FROM);
  boost::dynamic_bitset<> rev(width + HD_FROM);

  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto x: chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      int pos(chr.len -1 -x.F3);
      if(!RANGE(pos, start, end-1)) continue;
      if(strand==STRAND_PLUS) fwd.set(pos - start + HD_FROM);
      else                    rev.set(pos - start);
    }
  }
  for(int i=0; i<HD_WIDTH; ++i) {
    (fwd >>= 1);
    hd[i] += ((fwd ^ rev).count());
  }
  return;
}

void hammingDist(Mapfile &p, int numthreads)
{
  cout << "Making Hamming distance plot.." << flush;

#ifdef CHR1ONLY
  for(uint i=0; i<1; ++i) {
#else 
#pragma omp parallel for num_threads(numthreads)
  for(uint i=0; i<p.chr.size(); ++i) {
#endif
    hammingDistChr(p.chr[i], p.dist.hd);
  }

  //  GaussianSmoothing(p.dist.hd);

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

  long sum = accumulate(p.dist.hd.begin(), p.dist.hd.end(), 0.0);

  string filename = p.oprefix + ".hdp.csv";
  ofstream out(filename);
  out << "Strand shift\tHamming distance\tProp" << endl;
  for(int i=0; i<HD_WIDTH; ++i) out << (i - HD_FROM) << "\t" << p.dist.hd[i] << "\t" << (p.dist.hd[i]/(double)sum) << endl;

  cout << "done." << endl;
  cout << "Estimated fragment length: " << p.dist.eflen << endl;

  return;
}
  
