#include "pw_hamming.h"
#include "macro.h"
#include <map>
#include <boost/dynamic_bitset.hpp>
#include <boost/thread.hpp>
#include <math.h>   

//#define CHR1ONLY 1

enum functype {JACCARD, CCP, HAMMING};

void jaccardfunc_vector(Mapfile &p, SeqStats &chr, int numthreads);
void func_bitset(Mapfile &p, SeqStats &chr, functype type);
void ccpfunc(Mapfile &p, SeqStats &chr, int numthreads);

void shiftDist::funcCCP(vector<char> &fwd, vector<char> &rev, int start, int end)
{
  int width = end - start - ng_to;
  double mx,xx;
  double my,yy;
  calcMeanSD(fwd, width, mx, xx);
  calcMeanSD(rev, width, my, yy);
  
  for(int step=-mp_from; step<mp_to; step+=5) {
    double xy(0);
    for(int j=mp_from; j<width; ++j) xy += (fwd[j] - mx) * (rev[j+step] - my);
    mp[step] = xy;
  }
  
  for(int step=ng_from; step<ng_to; step+=ng_step) {
    double xy(0);
    for(int j=mp_from; j<width; ++j) xy += (fwd[j] - mx) * (rev[j+step] - my);
    nc[step] = xy;
  }

  double val = 1/(xx*yy*(width - mp_from - 1));
  for(auto itr = mp.begin(); itr != mp.end(); ++itr) itr->second *= val;
  for(auto itr = nc.begin(); itr != nc.end(); ++itr) itr->second *= val;
}

void shiftDist::funcJaccard(vector<char> &fwd, vector<char> &rev, int start, int end)
{
  int width = end - start - ng_to;
  int xx = accumulate(fwd.begin(), fwd.end(), 0);
  int yy = accumulate(rev.begin(), rev.end(), 0);
  
  for(int step=-mp_from; step<mp_to; step+=mp_step) {
    mp[step] = getJaccard(step, width, xx+yy, fwd, rev);
  }
  
  for(int step=ng_from; step<ng_to; step+=ng_step) {
    nc[step] = getJaccard(step, width, xx+yy, fwd, rev);
  }
}

void shiftDist::funcJaccard(boost::dynamic_bitset<> &fwd, const boost::dynamic_bitset<> &rev, int start, int end)
{
  int xx = fwd.count();
  int yy = rev.count();
  int xysum(xx+yy);
  
  fwd <<= mp_from;
  for(int step=-mp_from; step<mp_to; step+=mp_step) {
    fwd >>= 1;
    int xy((fwd & rev).count());
    mp[step] = xy/(double)(xysum-xy);
  }
  
  for(int step=ng_from; step<ng_to; ++step) {
    fwd >>= 1;
    int xy((fwd & rev).count());
    nc[step] = xy/(double)(xysum-xy);
  }
}

void shiftDist::funcHamming(boost::dynamic_bitset<> &fwd, const boost::dynamic_bitset<> &rev, int start, int end)
{
  fwd <<= mp_from;
  for(int step=-mp_from; step<mp_to; step+=mp_step) {
    fwd >>= 1;
    mp[step] = (fwd ^ rev).count();
  }
  
  for(int step=ng_from; step<ng_to; ++step) {
    fwd >>= 1;
    nc[step] = (fwd ^ rev).count();
  }
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

  boost::thread_group agroup;

#ifdef CHR1ONLY
  for(uint i=0; i<1; ++i) {
#else 
    //#pragma omp parallel for num_threads(numthreads)
  for(uint i=0; i<p.chr.size(); ++i) {
#endif
    agroup.create_thread(bind(jaccardfunc_vector, boost::ref(p), p.chr[i], numthreads));
    //    agroup.create_thread(bind(func_bitset, boost::ref(p), p.chr[i], JACCARD));
    //    jaccardfunc_vector(p, p.chr[i], numthreads);
    //func_bitset(p, p.chr[i], JACCARD);
  }

  agroup.join_all();

  string filename = p.oprefix + ".jaccard.csv";
  ofstream out(filename);
  outputmp(out, p.genome.jac, p, "Jaccard index");
  
  return;
}
  
void pw_ccp(Mapfile &p, int numthreads)
{
  printf("Making cross-correlation profile...\n");

  boost::thread_group agroup;
#ifdef CHR1ONLY
  for(uint i=0; i<1; ++i) {
#else 
    //#pragma omp parallel for num_threads(numthreads)
  for(uint i=0; i<p.chr.size(); ++i) {
#endif
    agroup.create_thread(bind(ccpfunc, boost::ref(p), p.chr[i], numthreads));
    //    ccpfunc(p, p.chr[i], numthreads);
  }
  
  string filename = p.oprefix + ".ccp.csv";
  ofstream out(filename);
  outputmp(out, p.genome.ccp, p, "Cross correlation");
  
  return;
}

void hammingDist(Mapfile &p, int numthreads)
{
  cout << "Making Hamming distance plot.." << flush;

  boost::thread_group agroup;
#ifdef CHR1ONLY
  for(uint i=0; i<1; ++i) {
#else 
    //#pragma omp parallel for num_threads(numthreads)
  for(uint i=0; i<p.chr.size(); ++i) {
#endif
    agroup.create_thread(bind(func_bitset, boost::ref(p), p.chr[i], HAMMING));
    //    func_bitset(p, p.chr[i], HAMMING);
  }

  //  GaussianSmoothing(p.dist.hd);

  /*  int min_hd_fl=p.dist.hd[HD_TO-1];
  int max_hd_fl=p.dist.hd[HD_TO-1];

  int bd = 2*p.dist.lenF3;
  for(int i=HD_TO-1; i>=bd+1; --i) {
    if(p.dist.hd[i] < min_hd_fl) {
      min_hd_fl = p.dist.hd[i];
      p.dist.eflen = i - HD_FROM;
    }
    if(max_hd_fl < p.dist.hd[i]) max_hd_fl = p.dist.hd[i];
    }*/

  string filename = p.oprefix + ".hdp.csv";
  ofstream out(filename);
  outputmp(out, p.genome.hd, p, "Hamming distance");

  return;
}

vector<char> genVector(const strandData &seq, int start, int end)
{
  vector<char> array((end-start), 0);
  for (auto x: seq.vRead) {
    if(!x.duplicate && RANGE(x.F3, start, end-1))
      ++array[x.F3 - start];
  }
  return array;
}

void ccpfunc(Mapfile &p, SeqStats &chr, int numthreads)
{
  cout << chr.name << endl;
  int start(0);
  int end(chr.len);
  //int start(1.5*NUM_100M);
  //int end(2*NUM_100M);
  //    int start(1.213*NUM_100M);
  //int end(1.214*NUM_100M);
  
  auto fwd = genVector(chr.seq[STRAND_PLUS], start, end);
  auto rev = genVector(chr.seq[STRAND_MINUS], start, end);

  chr.ccp.funcCCP(fwd, rev, start, end);

  if(chr.isautosome()) p.genome.addccp(chr);
  
  string filename = p.oprefix + ".ccp." + chr.name +".csv";
  ofstream out(filename);
  outputmp(out, chr.ccp, p, "Cross correlation");
  return;
}

double getJaccard(int step, int end, int xysum, const vector<char> &fwd, const vector<char> &rev)
{
  int xy(0);
  for(int j=HD_FROM; j<end; ++j) if(fwd[j] * rev[j+step]) xy += max(fwd[j], rev[j+step]);
  return (xy/(double)(xysum-xy));
  //    cout << xy << "\t" << xx << "\t" << yy << "\t" << endl;
  //  return;
}

void jaccardfunc_vector(Mapfile &p, SeqStats &chr, int numthreads)
{
  cout << chr.name << endl;
  //  chr.jac.setSE(1.213*NUM_100M, 1.214*NUM_100M);
  int start(0);
  int end(chr.len);
  
  auto fwd = genVector(chr.seq[STRAND_PLUS], start, end);
  auto rev = genVector(chr.seq[STRAND_MINUS], start, end);

  chr.jac.funcJaccard(fwd, rev, start, end);

  if(chr.isautosome()) p.genome.addjac(chr);
  
  string filename = p.oprefix + ".jaccard." + chr.name +".csv";
  ofstream out(filename);
  outputmp(out, chr.jac, p, "Jaccard index");
    
  return;
}

 boost::dynamic_bitset<> genBitset(const strandData &seq, int start, int end, long chrlen)
{
  boost::dynamic_bitset<> array(end - start + HD_FROM);
  for (auto x: seq.vRead) {
    if(x.duplicate) continue;
    int pos(chrlen -1 -x.F3);
    if(!RANGE(pos, start, end-1)) continue;
    array.set(pos - start);
  }
  return array;
}
 
void func_bitset(Mapfile &p, SeqStats &chr, functype type)
{
  cout << chr.name << endl;
  int start(0);
  int end(chr.len);
  //    int start(1.213*NUM_100M);
  //int end(1.214*NUM_100M);
  
  //  int width(end-start);
  
  auto fwd = genBitset(chr.seq[STRAND_PLUS], start, end, chr.len);
  auto rev = genBitset(chr.seq[STRAND_MINUS], start, end, chr.len);
  
  if(type==JACCARD) chr.jac.funcJaccard(fwd, rev, start, end);
  else if(type==HAMMING) chr.hd.funcHamming(fwd, rev, start, end);

  if(type==JACCARD) {
    if(chr.isautosome()) p.genome.addjac(chr);
    string filename = p.oprefix + ".jaccard." + chr.name +".csv";
    ofstream out(filename);
    outputmp(out, chr.jac, p, "Jaccard index");
  } else if(type==HAMMING) {
    if(chr.isautosome()) p.genome.addhd(chr);

    string filename = p.oprefix + ".hdp." + chr.name +".csv";
    ofstream out(filename);
    outputmp(out, chr.hd, p, "Hamming distance");
  }

  //  fwd <<= HD_FROM;
  /*  boost::dynamic_bitset<> fwd(width + HD_FROM);
  boost::dynamic_bitset<> rev(width + HD_FROM);

  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto x: chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      int pos(chr.len -1 -x.F3);
      if(!RANGE(pos, start, end-1)) continue;
      if(strand==STRAND_PLUS) fwd.set(pos - start + HD_FROM);
      else                    rev.set(pos - start);
    }
    }*/
  
  /*  int xx = fwd.count();
  int yy = rev.count();
  int xysum(xx+yy);

  for(int i=-HD_FROM; i<HD_TO; ++i) {
    (fwd >>= 1);
    if(type==JACCARD) {
      int xy((fwd & rev).count());
      chr.jac.mp[i] = xy/(double)(xysum-xy);
      //      cout << xy << "\t" << xx << "\t" << yy << "\t" << endl; 
    } else if(type==HAMMING) {
      chr.hd.mp[i] = ((fwd ^ rev).count());
    }
  }

  for(int i=NG_FROM; i<NG_TO; ++i) {
    (fwd >>= 1);
    if(type==JACCARD) {
      int xy((fwd & rev).count());
      chr.jac.nc[i] = xy/(double)(xysum-xy);
    } else if(type==HAMMING) {
      chr.hd.nc[i] = ((fwd ^ rev).count());
    }
    }*/

  return;
}
