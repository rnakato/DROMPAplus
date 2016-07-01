#include "pw_shiftprofile.h"
#include "pw_shiftprofile_p.h"
#include "macro.h"
#include <map>
#include <boost/dynamic_bitset.hpp>
#include <boost/thread.hpp>
#include <time.h>

//#define CHR1ONLY 1

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

void shiftCcp::setDist(int ichr, const vector<char> &fwd, const vector<char> &rev, int start, int end)
{
  int width = end - start - ng_to;
  double mx,xx;
  double my,yy;
  calcMeanSD(fwd, width, mx, xx);
  calcMeanSD(rev, width, my, yy);
  
  for(int step=-mp_from; step<mp_to; step+=5) {
    double xy(0);
    for(int j=mp_from; j<width; ++j) xy += (fwd[j] - mx) * (rev[j+step] - my);
    chr[ichr].mp[step] = xy;
  }
  
  for(int step=ng_from; step<ng_to; step+=ng_step) {
    double xy(0);
    for(int j=mp_from; j<width; ++j) xy += (fwd[j] - mx) * (rev[j+step] - my);
    chr[ichr].nc[step] = xy;
  }

  double val = 1/(xx*yy*(width - mp_from - 1));
  for(auto itr = chr[ichr].mp.begin(); itr != chr[ichr].mp.end(); ++itr) itr->second *= val;
  for(auto itr = chr[ichr].nc.begin(); itr != chr[ichr].nc.end(); ++itr) itr->second *= val;
}

double getJaccard(int step, int end, int xysum, const vector<char> &fwd, const vector<char> &rev)
{
  int xy(0);
  for(int j=HD_FROM; j<end; ++j) if(fwd[j] * rev[j+step]) xy += max(fwd[j], rev[j+step]);
  return (xy/(double)(xysum-xy));
}

void shiftJacVec::setDist(int ichr, const vector<char> &fwd, const vector<char> &rev, int start, int end)
{
  int width = end - start - ng_to;
  int xx = accumulate(fwd.begin(), fwd.end(), 0);
  int yy = accumulate(rev.begin(), rev.end(), 0);
  
  for(int step=-mp_from; step<mp_to; step+=mp_step) {
    chr[ichr].mp[step] = getJaccard(step, width, xx+yy, fwd, rev);
  }
  
  for(int step=ng_from; step<ng_to; step+=ng_step) {
    chr[ichr].nc[step] = getJaccard(step, width, xx+yy, fwd, rev);
  }
}

void shiftJacBit::setDist(int ichr, boost::dynamic_bitset<> &fwd, const boost::dynamic_bitset<> &rev, int start, int end)
{
  int xx = fwd.count();
  int yy = rev.count();
  int xysum(xx+yy);
  
  fwd <<= mp_from;
  for(int step=-mp_from; step<mp_to; step+=mp_step) {
    fwd >>= 1;
    int xy((fwd & rev).count());
    chr[ichr].mp[step] = xy/(double)(xysum-xy);
  }
  
  for(int step=ng_from; step<ng_to; ++step) {
    fwd >>= 1;
    int xy((fwd & rev).count());
    chr[ichr].nc[step] = xy/(double)(xysum-xy);
  }
}

void shiftHamming::setDist(int ichr, boost::dynamic_bitset<> &fwd, const boost::dynamic_bitset<> &rev, int start, int end)
{
  fwd <<= mp_from;
  for(int step=-mp_from; step<mp_to; step+=mp_step) {
    fwd >>= 1;
    chr[ichr].mp[step] = (fwd ^ rev).count();
  }
  
  for(int step=ng_from; step<ng_to; ++step) {
    fwd >>= 1;
    chr[ichr].nc[step] = (fwd ^ rev).count();
  }
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


template <class T>
void genThread(T &dist, const Mapfile &p, uint s, uint e, string typestr, boost::mutex &mtx) {
  for(uint i=s; i<=e; ++i) {
    cout << p.chr[i].name << endl;
    int start(0);
    int end(p.chr[i].len);
  //    int start(1.213*NUM_100M);
  //int end(1.214*NUM_100M);

    dist.execchr(p, i, start, end, mtx);

    if(p.chr[i].isautosome()) dist.add2genome(dist.chr[i], mtx);
 
    string filename = p.oprefix + "." + typestr + "." + p.chr[i].name + ".csv";
    dist.chr[i].outputmp(filename, dist.genome.nread, dist.name);
  }
}

template <class T>
void makeProfile(Mapfile &p, string typestr, int numthreads)
{
  T dist(p);
  cout << "Making " << dist.name << " profile..." << endl;

  boost::thread_group agroup;
  boost::mutex mtx;
  
#ifdef CHR1ONLY
  genThread(dist, p, 0, 0, typestr, mtx);
#else 
  for(uint i=0; i<p.vsepchr.size(); i++) {
    agroup.create_thread(bind(genThread<T>, boost::ref(dist), boost::cref(p), p.vsepchr[i].s, p.vsepchr[i].e, typestr, boost::ref(mtx)));
  }
  agroup.join_all();
#endif
  
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

  string filename = p.oprefix + "." + typestr + ".csv";
  dist.genome.outputmp(filename, dist.genome.nread, dist.name);
  
  return;
}

void strShiftProfile(Mapfile &p, string typestr, int numthreads)
{
  if(typestr=="exjaccard") makeProfile<shiftJacVec>(p, typestr, numthreads);
  if(typestr=="jaccard")   makeProfile<shiftJacBit>(p, typestr, numthreads);
  else if(typestr=="ccp")  makeProfile<shiftCcp>(p,    typestr, numthreads);
  else if(typestr=="hdp")  makeProfile<shiftHamming>(p, typestr, numthreads);
  
  return;
}
