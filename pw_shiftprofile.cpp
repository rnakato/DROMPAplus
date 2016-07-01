#include "pw_shiftprofile.h"
#include "pw_shiftprofile_p.h"
#include "macro.h"
#include <map>
#include <boost/dynamic_bitset.hpp>
#include <boost/thread.hpp>
#include <time.h>

//#define CHR1ONLY 1

double getJaccard(int step, int end, int xysum, int mp_from, const vector<char> &fwd, const vector<char> &rev)
{
  int xy(0);
  for(int j=mp_from; j<end; ++j) if(fwd[j] * rev[j+step]) xy += max(fwd[j], rev[j+step]);
  return (xy/(double)(xysum-xy));
}

void genThreadJacVec(_shiftDist &chr, int xysum, int mp_from, const vector<char> &fwd, const vector<char> &rev, int s, int e, boost::mutex &mtx)
{
  for(int step=s; step<e; ++step) {
    chr.setmp(step, getJaccard(step, chr.width4mp, xysum, mp_from, fwd, rev), mtx);
  }
}

void shiftJacVec::setDist(_shiftDist &chr, const vector<char> &fwd, const vector<char> &rev)
{
  int xx = accumulate(fwd.begin(), fwd.end(), 0);
  int yy = accumulate(rev.begin(), rev.end(), 0);

  boost::thread_group agroup;
  boost::mutex mtx;
  for(uint i=0; i<seprange.size(); i++) {
    agroup.create_thread(bind(&genThreadJacVec, boost::ref(chr), xx+yy, mp_from, boost::cref(fwd), boost::cref(rev), seprange[i].start, seprange[i].end, boost::ref(mtx)));
  }
  agroup.join_all();

  for(int step=ng_from; step<ng_to; step+=ng_step) {
    chr.nc[step] = getJaccard(step, chr.width4mp, xx+yy, mp_from, fwd, rev);
  }
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

void genThreadCcp(_shiftDist &chr, const vector<char> &fwd, const vector<char> &rev, double mx, double my, int mp_from, const int s, const int e, boost::mutex &mtx)
{
  for(int step=s; step<e; ++step) {
    double xy(0);
    for(int j=mp_from; j<chr.width4mp; ++j) xy += (fwd[j] - mx) * (rev[j+step] - my);
    chr.setmp(step, xy, mtx);
  }
}

void shiftCcp::setDist(_shiftDist &chr, const vector<char> &fwd, const vector<char> &rev)
{
  double mx, my, xx, yy;
  calcMeanSD(fwd, chr.width4mp, mx, xx);
  calcMeanSD(rev, chr.width4mp, my, yy);

  boost::thread_group agroup;
  boost::mutex mtx;
  for(uint i=0; i<seprange.size(); i++) {
    agroup.create_thread(bind(genThreadCcp, boost::ref(chr), boost::cref(fwd), boost::cref(rev), mx, my, mp_from, seprange[i].start, seprange[i].end, boost::ref(mtx)));
  }
  agroup.join_all();
  
  for(int step=ng_from; step<ng_to; step+=ng_step) {
    double xy(0);
    for(int j=mp_from; j<chr.width4mp; ++j) xy += (fwd[j] - mx) * (rev[j+step] - my);
    chr.nc[step] = xy;
  }

  double val = 1/(xx*yy*(chr.width4mp - mp_from - 1));
  for(auto itr = chr.mp.begin(); itr != chr.mp.end(); ++itr) itr->second *= val;
  for(auto itr = chr.nc.begin(); itr != chr.nc.end(); ++itr) itr->second *= val;
}

void shiftJacBit::setDist(_shiftDist &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev)
{
  int xysum(fwd.count() + rev.count());

  rev <<= mp_from;
  for(int step=-mp_from; step<mp_to; ++step) {
    rev >>= 1;
    int xy((fwd & rev).count());
    chr.mp[step] = xy/(double)(xysum-xy);
  }
  rev >>= (ng_from - mp_to);
  
  for(int step=ng_from; step<ng_to; step+=ng_step) {
    rev >>= ng_step;
    int xy((fwd & rev).count());
    chr.nc[step] = xy/(double)(xysum-xy);
  }
}

void shiftHamming::setDist(_shiftDist &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev)
{
  rev <<= mp_from;
  for(int step=-mp_from; step<mp_to; ++step) {
    rev >>= 1;
    chr.mp[step] = (fwd ^ rev).count();
  }
  rev >>= (ng_from - mp_to);
  
  for(int step=ng_from; step<ng_to; step+=ng_step) {
    rev >>= ng_step;
    chr.nc[step] = (fwd ^ rev).count();
  }
}

vector<char> genVector(const strandData &seq, int start, int end)
{
  vector<char> array(end-start, 0);
  for (auto x: seq.vRead) {
    if(!x.duplicate && RANGE(x.F3, start, end-1))
      ++array[x.F3 - start];
  }
  return array;
}

boost::dynamic_bitset<> genBitset(const strandData &seq, int start, int end)
{
  boost::dynamic_bitset<> array(end-start);
  for (auto x: seq.vRead) {
    if(!x.duplicate && RANGE(x.F3, start, end-1))
      array.set(x.F3 - start);
  }
  return array;
}

template <class T>
void genThread(T &dist, const Mapfile &p, uint s, uint e, string typestr, boost::mutex &mtx) {
  for(uint i=s; i<=e; ++i) {
    cout << p.chr[i].name << endl;
   
    dist.execchr(p, i);

    if(p.chr[i].isautosome()) dist.add2genome(dist.chr[i], mtx);
 
    string filename = p.oprefix + "." + typestr + "." + p.chr[i].name + ".csv";
    dist.chr[i].outputmp(filename, dist.genome.nread, dist.name);
  }
}

template <class T>
void makeProfile(Mapfile &p, string typestr, int numthreads)
{
  T dist(p, numthreads);
  cout << "Making " << dist.name << " profile..." << endl;

  boost::thread_group agroup;
  boost::mutex mtx;

#ifdef CHR1ONLY
  genThread(dist, p, 0, 0, typestr, mtx);
#else 
  if(typestr == "hdp" || typestr == "jaccard") {
    for(uint i=0; i<p.vsepchr.size(); i++) {
      agroup.create_thread(bind(genThread<T>, boost::ref(dist), boost::cref(p), p.vsepchr[i].s, p.vsepchr[i].e, typestr, boost::ref(mtx)));
    }
    agroup.join_all();
  } else {
    genThread(dist, p, 0, p.chr.size()-1, typestr, mtx);
  }
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
