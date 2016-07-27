#include "pw_shiftprofile.h"
#include "pw_shiftprofile_p.h"
#include "macro.h"
#include <map>
#include <boost/dynamic_bitset.hpp>
#include <boost/thread.hpp>

/*namespace {
  enum {BP_BACKGROUD, BP_REPEAT, BP_PEAK};
  }*/

/*void ReadShiftProfileAll::defSepRange(int numthreads)
{
  int length(mp_to+mp_from);
  int sepsize = length/numthreads +1;
  for(int i=0; i<numthreads; i++) {
    int s = i*sepsize;
    int e = (i+1)*sepsize;
    if(i==numthreads-1) e = length;
    range sep(s - mp_from, e - mp_from);
    seprange.push_back(sep);
  }
  }*/

void addmp(std::map<int, double> &mpto, const std::map<int, double> &mpfrom, double w)
{
  for(auto itr = mpfrom.begin(); itr != mpfrom.end(); ++itr) {
    mpto[itr->first] += itr->second * w;
  }
}

double getJaccard(int step, int width, int xysum, const std::vector<char> &fwd, const std::vector<char> &rev)
{
  int xy(0);
  for(int j=mp_from; j<width-ng_to; ++j) if(fwd[j] * rev[j+step]) xy += std::max(fwd[j], rev[j+step]);
  return (xy/static_cast<double>(xysum-xy));
}

void genThreadJacVec(ReadShiftProfile &chr, int xysum, const std::vector<char> &fwd, const std::vector<char> &rev, int s, int e, boost::mutex &mtx)
{
  for(int step=s; step<e; ++step) {
    chr.setmp(step, getJaccard(step, chr.width, xysum, fwd, rev), mtx);
  }
}

void shiftJacVec::setDist(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev)
{
  int xx = accumulate(fwd.begin(), fwd.end(), 0);
  int yy = accumulate(rev.begin(), rev.end(), 0);

  boost::thread_group agroup;
  boost::mutex mtx;
  for(uint i=0; i<seprange.size(); i++) {
    agroup.create_thread(bind(&genThreadJacVec, boost::ref(chr), xx+yy, boost::cref(fwd), boost::cref(rev), seprange[i].start, seprange[i].end, boost::ref(mtx)));
  }
  agroup.join_all();

  for(int step=ng_from; step<ng_to; step+=ng_step) {
    chr.nc[step] = getJaccard(step, chr.width-ng_to, xx+yy, fwd, rev);
  }
}

void genThreadCcp(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev, double mx, double my, const int s, const int e, boost::mutex &mtx)
{
  for(int step=s; step<e; ++step) {
    double xy(0);
    for(int j=mp_from; j<chr.width-ng_to; ++j) xy += (fwd[j] - mx) * (rev[j+step] - my);
    chr.setmp(step, xy, mtx);
  }
}

void shiftCcp::setDist(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev)
{
  moment<char> x(fwd, mp_from, chr.width-ng_to);
  moment<char> y(rev, mp_from, chr.width-ng_to);

  boost::thread_group agroup;
  boost::mutex mtx;
  for(uint i=0; i<seprange.size(); i++) {
    agroup.create_thread(bind(genThreadCcp, boost::ref(chr), boost::cref(fwd), boost::cref(rev), x.getmean(), y.getmean(), seprange[i].start, seprange[i].end, boost::ref(mtx)));
  }
  agroup.join_all();

  for(int step=ng_from; step<ng_to; step+=ng_step) {
    double xy(0);
    for(int j=mp_from; j<chr.width-ng_to; ++j) xy += (fwd[j] - x.getmean()) * (rev[j+step] - y.getmean());
    chr.nc[step] = xy;
  }

  double val = 1/(x.getsd() * y.getsd() * (chr.width-ng_to - mp_from - 1));
  for(auto itr = chr.mp.begin(); itr != chr.mp.end(); ++itr) itr->second *= val;
  for(auto itr = chr.nc.begin(); itr != chr.nc.end(); ++itr) itr->second *= val;
}

void shiftJacBit::setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev)
{
  int xysum(fwd.count() + rev.count());

  rev <<= mp_from;
  for(int step=-mp_from; step<mp_to; ++step) {
    rev >>= 1;
    int xy((fwd & rev).count());
    chr.mp[step] = xy/static_cast<double>(xysum-xy);
  }
  rev >>= (ng_from - mp_to);
  
  for(int step=ng_from; step<ng_to; step+=ng_step) {
    rev >>= ng_step;
    int xy((fwd & rev).count());
    chr.nc[step] = xy/static_cast<double>(xysum-xy);
  }
}

void shiftHamming::setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev)
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

std::vector<char> genVector(const strandData &seq, int start, int end)
{
  std::vector<char> array(end-start, 0);
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
void genThread(T &dist, const Mapfile &p, uint chr_s, uint chr_e, std::string typestr) {
  for(uint i=chr_s; i<=chr_e; ++i) {
    std::cout << p.chr[i].name << ".." << std::flush;
   
    dist.execchr(p, i);
    dist.chr[i].setflen();
  }
}

template <class T>
int getRepeatRegion(std::vector<range> &vrep, int j, std::vector<T> array, int start, int end)
{
  int lower_thre=5;
  int s, e;
  for(s=j; s>=start; --s) if(array[s]<lower_thre) break;
  for(e=j; e<end; ++e)    if(array[e]<lower_thre) break;

  vrep.push_back(range(s,e));
  
  return e;
}

/*template <class T>
void func(T &dist, const Mapfile &p, const int i) {

  auto fwd = genBitset(p.chr[i].seq[STRAND_PLUS],  dist.chr[i].start, dist.chr[i].end);
  auto rev = genBitset(p.chr[i].seq[STRAND_MINUS], dist.chr[i].start, dist.chr[i].end);

  int flen(dist.getnsci());
  int readlen(p.getlenF3());
  
  std::vector<short> fragarray(dist.chr[i].width, 0);
  std::vector<short> reparray(dist.chr[i].width, 0);
  for(int j=dist.chr[i].start; j<dist.chr[i].end-flen; ++j) {
    if(fwd[j] && rev[j+flen]) {
      ++dist.chr[i].numOfFragmentWithFlen;
      for(int k=0; k<flen; ++k) ++fragarray[j - dist.chr[i].start +k];
    }
    if(fwd[j] && rev[j+readlen]) {
      ++dist.chr[i].numOfFragmentWithReplen;
      for(int k=0; k<readlen; ++k) ++reparray[j - dist.chr[i].start +k];
    }
  }

  for(int j=dist.chr[i].start; j<dist.chr[i].end-flen; ++j) {
    if(fragarray[j]) ++dist.chr[i].numOfCoveredBaseWithFlen;
    if(reparray[j])  ++dist.chr[i].numOfCoveredBaseWithReplen;
  }

  std::vector<int> dfrag(10000,0);
  std::vector<int> drep(10000,0);
  for(int j=dist.chr[i].start; j<dist.chr[i].end; ++j) {
    ++dfrag[fragarray[j]];
    ++drep[reparray[j]];
  }

  double pdfrag(0), pdrep(0);
  double dfragsum = accumulate(dfrag.begin(), dfrag.end(), 0);
  double drepsum = accumulate(drep.begin(), drep.end(), 0);
  int thre4fragarray(0);
  for(int j=0; j<flen; ++j) {
    double a(dfrag[j]/dfragsum);
    pdfrag += a;
    if(!thre4fragarray && pdfrag > 0.95) thre4fragarray = j;
    double b(drep[j]/drepsum);
    pdrep += b;
    std::cerr << j << "\t" << a << "\t" << pdfrag << "\t" << b<< "\t" << pdrep << std::endl;
  }

  std::cout << flen << " thre4fragarray " << thre4fragarray << std::endl;
  
  std::vector<range> vrep;
  int thre4reparray(10);
  for(int j=dist.chr[i].start; j<dist.chr[i].end; ++j) {
    if(reparray[j]>=thre4reparray) j = getRepeatRegion(vrep, j, reparray, dist.chr[i].start, dist.chr[i].end);
  }

  std::vector<char> anoarray(dist.chr[i].width, BP_BACKGROUD);
  for(int j=dist.chr[i].start; j<dist.chr[i].end; ++j) {
    if(fragarray[j] > thre4fragarray) anoarray[j-dist.chr[i].start] = BP_PEAK;
  }
  for(auto x:vrep) {
    for(int j=x.start; j<x.end; ++j) anoarray[j] = BP_REPEAT;
  }

  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto x: p.chr[i].seq[strand].vRead) {
      if(x.duplicate || !RANGE(x.F3, dist.chr[i].start, dist.chr[i].end-1)) continue;
      if(anoarray[x.F3 - dist.chr[i].start]==BP_BACKGROUD)   ++dist.chr[i].nread_back;
      else if(anoarray[x.F3 - dist.chr[i].start]==BP_REPEAT) ++dist.chr[i].nread_rep;
      else ++dist.chr[i].nread_peak;
    }
  }

  //  std::cout << dist.chr[i].nread << "\t" << dist.chr[i].nread_back << "\t" << dist.chr[i].nread_peak << "\t" << dist.chr[i].nread_rep << std::endl;

  return;
}
*/

template <class T>
void func(T &dist, const Mapfile &p, const int i) {
  auto fwd = genBitset(p.chr[i].seq[STRAND_PLUS],  dist.chr[i].start, dist.chr[i].end);
  auto rev = genBitset(p.chr[i].seq[STRAND_MINUS], dist.chr[i].start, dist.chr[i].end);

  int flen(dist.getnsci());
  int readlen(p.getlenF3());

  dist.chr[i].setFragmentVariability4Frag(flen, fwd, rev);
  dist.chr[i].setFragmentVariability4Rep(readlen, fwd, rev);
  dist.chr[i].setFragmentVariability4Back(ng_to, fwd, rev);

  return;
}

template <class T>
void genThread_countbkreads(T &dist, const Mapfile &p, uint chr_s, uint chr_e, std::string typestr, boost::mutex &mtx)
{
  for(uint i=chr_s; i<=chr_e; ++i) {
    std::cout << p.chr[i].name << ".." << std::flush;
    func(dist, p, i);
    if(p.chr[i].isautosome()) dist.addnread2genome(i, mtx);
 
    std::string filename = p.getprefix() + "." + typestr + "." + p.chr[i].name + ".csv";
    dist.outputmpChr(filename, i);
  }
}

template <class T>
void makeProfile(Mapfile &p, const std::string &typestr, const int numthreads)
{
  T dist(p, numthreads);
  dist.printStartMessage();

  boost::thread_group agroup;
  boost::mutex mtx;

  if(typestr == "hdp" || typestr == "jaccard") {
    for(uint i=0; i<p.vsepchr.size(); i++) {
      agroup.create_thread(bind(genThread<T>, boost::ref(dist), boost::cref(p), p.vsepchr[i].s, p.vsepchr[i].e, typestr));
    }
    agroup.join_all();
  } else {
    //genThread(dist, p, 0, p.chr.size()-1, typestr);
    genThread(dist, p, 0, 0, typestr);
  }

  //  GaussianSmoothing(p.dist.hd);

  // set fragment length;
  for(uint i=0; i<p.chr.size(); ++i) {
    if(p.chr[i].isautosome()) dist.addmp2genome(i);
  }
  

  dist.setflen();
  p.seteflen(dist.getnsci());

  std::cout << "count reads in background.." << std::flush;
  for(uint i=0; i<p.vsepchr.size(); i++) {
    agroup.create_thread(bind(genThread_countbkreads<T>, boost::ref(dist), boost::cref(p), p.vsepchr[i].s, p.vsepchr[i].e, typestr, boost::ref(mtx)));
  }
  agroup.join_all();
  
  std::string filename = p.getprefix() + "." + typestr + ".csv";
  dist.outputmpGenome(filename);

  return;
}

void strShiftProfile(Mapfile &p, std::string typestr, int numthreads)
{
  if(typestr=="exjaccard") makeProfile<shiftJacVec>(p, typestr, numthreads);
  if(typestr=="jaccard")   makeProfile<shiftJacBit>(p, typestr, numthreads);
  else if(typestr=="ccp")  makeProfile<shiftCcp>(p, typestr, numthreads);
  else if(typestr=="hdp")  makeProfile<shiftHamming>(p, typestr, numthreads);
  
  return;
}
