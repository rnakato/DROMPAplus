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
  //  cout << chr.name <<".." << flush;
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
  
#pragma omp parallel for num_threads(numthreads)
#ifdef CHR1ONLY
  for(uint i=0; i<1; ++i) {
#else 
  for(uint i=0; i<p.chr.size(); ++i) {
#endif
    hammingDistChr(p.chr[i], p.dist.hd, numthreads);
  }
  /*  if(p.lchr->len > NUM_10M) hammingDistChr(*p.lchr, p.dist.hd, numthreads);
  else {
    for(auto &chr: p.chr) hammingDistChr(chr, p.dist.hd, numthreads);
    }*/

  //  GaussianSmoothing(p.dist.hd);

   // get fragment length FL and HD[FL] run through from (i_num-1),...,2*read_len+1
  int min_hd_fl=p.dist.hd[HD_WIDTH-1];
  int max_hd_fl=p.dist.hd[HD_WIDTH-1];

  //  cout << p.dist.lenF3 << endl;
  int bd = 2*p.dist.lenF3;
  for(int i=HD_WIDTH-1; i>=bd+1; --i) {
    if(p.dist.hd[i] < min_hd_fl) {
      min_hd_fl = p.dist.hd[i];
      p.dist.eflen = i - HD_FROM;
    }
    if(max_hd_fl < p.dist.hd[i]) max_hd_fl = p.dist.hd[i];
  }

  //  long sum = accumulate(p.dist.hd.begin(), p.dist.hd.end(), 0);
  long sum(0);
  for(uint i=0; i<p.dist.hd.size(); ++i) sum += p.dist.hd[i];

  string filename = p.oprefix + ".hdp.csv";
  ofstream out(filename);
  out << "Strand shift\tHamming distance\tProp" << endl;
  for(int i=0; i<HD_WIDTH; ++i) out << (i - HD_FROM) << "\t" << p.dist.hd[i] << "\t" << (p.dist.hd[i]/(double)sum) << endl;

  cout << "done." << endl;
  cout << "Estimated fragment length: " << p.dist.eflen << endl;

  return;
}


double getControlRatio(map<int, double> &mp)
{
  double r(0);
  int n(0);
  for(auto itr = mp.begin(); itr != mp.end(); ++itr) {
    r += itr->second; 
    ++n;
    cout << itr->second << endl;
  }
  r /= n;
  r = 1/r;
  cout << r << endl;
  return r;
}
 
 void outputMap(ofstream &out, map<int, double> &mp, map<int, double> &mpnc, string str)
{
  double sum(0);
  for(auto itr = mp.begin(); itr != mp.end(); ++itr) sum += itr->second;

  double r = getControlRatio(mpnc);
  double nsc(0);
  int nsci(0);
  for(auto itr = mp.begin(); itr != mp.end(); ++itr) {
    if(nsc < itr->second*r) {
      nsc = itr->second*r;
      nsci = itr->first;
    }
  }
  
  out << "NSC\t"<< nsc << endl;
  out << str << "\t" << nsci << endl;
  out << "Strand shift\tCross correlation\tprop\tper control" << endl;
  for(auto itr = mp.begin(); itr != mp.end(); ++itr) {
    out << itr->first << "\t" << itr->second << "\t" << (itr->second/sum)<< "\t" << (itr->second*r) << endl;
  }
  return;
}

 double getJaccard(vector<char> fwd, vector<char> rev, int step, int xx, int yy, int max, double r, int numthreads)
{
  int xy(0);
#pragma omp parallel for num_threads(numthreads)
  for(int j=HD_FROM; j<max; ++j) {
    xy += fwd[j] * rev[j+step];
    //    cout << step<< "\t" << (xy/(double)(xx+yy-xy)) << "\t" << xy << "\t" << xx<< "\t" << yy << endl;
  }
  double val = xy/(double)(xx+yy-xy) * r;
  return val;
}
 
void pw_Jaccard(Mapfile &p, int numthreads)
{
  printf("Making Jaccard index profile...\n");

  map<int, double> mp;
  for(int i=-HD_FROM; i<HD_WIDTH; i+=5) mp[i] = 0;
  map<int, double> mpnc;
  for(int i=HD_NG_FROM; i<HD_NG_TO; i+=HD_NG_STEP) mpnc[i] = 0;

  for (auto chr: p.chr) {
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

    int max = width - HD_WIDTH;
    int xx(0), yy(0);
    for(int i=HD_FROM; i<max; ++i) {
      xx += fwd[i];
      yy += rev[i];
    }

#pragma omp parallel for num_threads(numthreads)
    for(int step=-HD_FROM; step<HD_WIDTH; ++step) {
      mp[step] += getJaccard(fwd, rev, step, xx, yy, max, chr.bothnread_nonred()/(double)p.genome.bothnread_nonred(), numthreads);
    }
    
    for(int step=HD_NG_FROM; step<HD_NG_TO; step+=HD_NG_STEP) {
      mpnc[step] += getJaccard(fwd, rev, step, xx, yy, max, chr.bothnread_nonred()/(double)p.genome.bothnread_nonred(), numthreads);
    }

#ifdef CHR1ONLY
    break;
#endif
  }

  string filename = p.oprefix + ".jaccard.csv";
  ofstream out(filename);
  outputMap(out, mp, mpnc, "Jaccard index");
  
  return;
}

void calcMeanSD(char *x, int max, double &ave, double &sd)
{
  double dx, var(0);

  ave=0;
  for(int i=HD_FROM; i<max; ++i) ave += x[i];
  ave /= (double)max;
  for(int i=HD_FROM; i<max; ++i) {
    dx = x[i] - ave;
    var += dx * dx;
  }
  sd = sqrt(var/double(max-1));
}
 
void pw_ccp(Mapfile &p, int numthreads)
{
  printf("Making cross-correlation profile...\n");

  map<int, double> mp;
  for(int i=-HD_FROM; i<HD_WIDTH; i+=5) mp[i] = 0;
  map<int, double> mpnc;
  for(int i=HD_NG_FROM; i<HD_NG_TO; i+=HD_NG_STEP) mpnc[i] = 0;
  
  for (auto chr: p.chr) {
    cout << chr.name << endl;
    int start(0);
    int end(chr.len);
    //int start(1.5*NUM_100M);
    //int end(2*NUM_100M);
    //    int start(1.213*NUM_100M);
    //int end(1.214*NUM_100M);
    
    int width(end-start);

    char *fwd = (char *)calloc(width, sizeof(char));
    char *rev = (char *)calloc(width, sizeof(char));

    for(int strand=0; strand<STRANDNUM; ++strand) {
      for (auto x: chr.seq[strand].vRead) {
	if(x.duplicate) continue;
	if(!RANGE(x.F3, start, end-1)) continue;
	if(strand==STRAND_PLUS) ++fwd[x.F3 - start];
	else                    ++rev[x.F3 - start];
      }
    }

    int max = width - HD_WIDTH;
    double mx,xx;
    double my,yy;
    calcMeanSD(fwd, max, mx, xx);
    calcMeanSD(rev, max, my, yy);

#pragma omp parallel for num_threads(numthreads)
    for(int step=-HD_FROM; step<HD_WIDTH; step+=5) {
      double xy(0);
#pragma omp parallel for num_threads(numthreads) reduction(+:xy)
      for(int j=HD_FROM; j<max; ++j) {
	xy += (fwd[j] - mx) * (rev[j+step] - my);
      }
      xy /= (double)(max - HD_FROM - 1);
      mp[step] += (xy / (xx*yy)) * (chr.bothnread_nonred()/(double)p.genome.bothnread_nonred());
    }

    for(int step=HD_NG_FROM; step<HD_NG_TO; step+=HD_NG_STEP) {
      double xy(0);
#pragma omp parallel for num_threads(numthreads) reduction(+:xy)
      for(int j=HD_FROM; j<max; ++j) {
	xy += (fwd[j] - mx) * (rev[j+step] - my);
      }
      xy /= (double)(max - HD_FROM - 1);
      mpnc[step] += (xy / (xx*yy)) * (chr.bothnread_nonred()/(double)p.genome.bothnread_nonred());
    }
    
    free(fwd);
    free(rev);
#ifdef CHR1ONLY
    break;
#endif
  }

  string filename = p.oprefix + ".ccp.csv";
  ofstream out(filename);
  outputMap(out, mp, mpnc, "Estimated fragment length");
  
  return;
}

