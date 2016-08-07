#include "pw_shiftprofile.h"
#include "pw_shiftprofile_p.h"
#include "macro.h"
#include <map>
#include <boost/thread.hpp>

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
    //    chr.mp[step] = xy;
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

std::vector<char> genVector4FixedReadsNum(const strandData &seq, int start, int end, const double r4cmp)
{
  static int n(0);
  int nseq(0);
  std::vector<char> array(end-start, 0);
  for (auto x: seq.vRead) {
    if(!x.duplicate && RANGE(x.F3, start, end-1)){
    //    if(!x.duplicate && RANGE(x.F3, start, 120000000)) {
      if(rand() >= r4cmp) continue;
      ++array[x.F3 - start];
      //      ++n;
      //++nseq;
    }
  }

  //  std::cout << "total num: " << n  << ", each: " << nseq << std::endl;
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
    std::cout << p.genome.chr[i].name << ".." << std::flush;

    dist.execchr(p, i);
    dist.chr[i].setflen();
    
    std::string filename = p.getprefix() + "." + typestr + "." + p.genome.chr[i].name + ".csv";
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
    //agroup.create_thread(bind(genThread<T>, boost::ref(dist), boost::cref(p), 0, 0, typestr));
    for(uint i=0; i<p.genome.vsepchr.size(); i++) {
      agroup.create_thread(bind(genThread<T>, boost::ref(dist), boost::cref(p), p.genome.vsepchr[i].s, p.genome.vsepchr[i].e, typestr));
    }
    agroup.join_all();
  } else {
    genThread(dist, p, 0, p.genome.chr.size()-1, typestr);
    //genThread(dist, p, 0, 0, typestr);
  }

  // set fragment length;
  for(uint i=0; i<p.genome.chr.size(); ++i) {
    if(p.genome.chr[i].isautosome()) dist.addmp2genome(i);
  }

  dist.setflen();
  p.seteflen(dist.getnsci());

  std::string filename = p.getprefix() + "." + typestr + ".csv";
  dist.outputmpGenome(filename);

  return;
}

void makeFragVarProfile(Mapfile &p, const std::string &typestr, const int numthreads, const int numRead4fvp, const int flen, const bool fvpfull)
{
  shiftFragVar dist(p, numthreads, flen, fvpfull);
  dist.printStartMessage();
  
  double r(1);
  if(numRead4fvp) r = numRead4fvp/static_cast<double>(dist.getnread());
  if(r>1){
    std::cerr << "\nWarning: number of reads for Fragment variability is < "<< (int)(numRead4fvp/NUM_1M) <<" million.\n";
    dist.lackOfReads_on();
  }

  double r4cmp = r*RAND_MAX;

  for(uint i=0; i<=p.genome.chr.size()-1; ++i) {
    if(p.genome.chr[i].isautosome()) {
      std::cout << p.genome.chr[i].name << ".." << std::flush;
      dist.execchr(p, i, r4cmp);
      std::string filename = p.getprefix() + "." + typestr + "." + p.genome.chr[i].name + ".csv";
      dist.outputmpChr(filename, i);
      dist.addmp2genome(i);
    }
  }

  std::string filename1 = p.getprefix() + ".mpfv.csv";
  dist.printmpfv(filename1);
  std::string filename2 = p.getprefix() + "." + typestr + ".csv";
  dist.outputmpGenome(filename2);

  return;
}

void strShiftProfile(const MyOpt::Variables &values, Mapfile &p, std::string typestr)
{
  int numthreads(values["threads"].as<int>());
  if(typestr=="exjaccard")    makeProfile<shiftJacVec>(p, typestr, numthreads);
  else if(typestr=="jaccard") makeProfile<shiftJacBit>(p, typestr, numthreads);
  else if(typestr=="ccp")     makeProfile<shiftCcp>(p, typestr, numthreads);
  else if(typestr=="hdp")     makeProfile<shiftHamming>(p, typestr, numthreads);
  else if(typestr=="fvp")     makeFragVarProfile(p, typestr, numthreads, values["nfvp"].as<int>(), p.getflen(values), values.count("fvpfull"));
  
  return;
}

void genThreadFragVar(ReadShiftProfile &chr, std::map<int, FragmentVariability> &mpfv, const std::vector<char> &fwd, const std::vector<char> &rev, const std::vector<double> &fvback, const int s, const int e, boost::mutex &mtx)
{
  for(int step=s; step<e; ++step) {
    FragmentVariability fv;
    fv.setVariability(step, chr.start, chr.end, fwd, rev);

    double diffMax(0);
    for(size_t k=0; k<sizeOfvDistOfDistaneOfFrag; ++k) {
      //      std::cout << fv.getAccuOfDistanceOfFragment(k) << "\t" << fvback[k] << std::endl;
      diffMax = std::max(diffMax, fv.getAccuOfDistanceOfFragment(k) - fvback[k]);
    }
    chr.setmp(step, diffMax, mtx);
    mpfv[step].add2genome(fv, mtx);
  }
}

void shiftFragVar::setDist(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev)
{
  boost::thread_group agroup;
  boost::mutex mtx;

  std::vector<double> fvback(sizeOfvDistOfDistaneOfFrag,0);
  int n(0);
  //  for(int step=ng_from; step<ng_to; step+=ng_step) {
  for(int step=NUM_100K; step<NUM_1M; step+=NUM_100K) {
    FragmentVariability fv;
    fv.setVariability(step, chr.start, chr.end, fwd, rev);    
    for(size_t k=0; k<sizeOfvDistOfDistaneOfFrag; ++k) {
      fvback[k] += fv.getAccuOfDistanceOfFragment(k);
    }
    ++n;

  }
  for(size_t k=0; k<sizeOfvDistOfDistaneOfFrag; ++k) fvback[k] /= n;

  if (fvpfull) {
    for(uint i=0; i<seprange.size(); i++) {
      agroup.create_thread(bind(&genThreadFragVar, boost::ref(chr), boost::ref(mpfv), boost::cref(fwd), boost::cref(rev), boost::cref(fvback), seprange[i].start, seprange[i].end, boost::ref(mtx)));
    }
    agroup.join_all();
  } else {
    //    std::vector<int> v{flen, chr.getlenF3()};
    // std::copy(v4mpfv.begin(), v4mpfv.end(), std::back_inserter(v));
    // for(auto x: v) {
    for(auto x: v4mpfv) {
      FragmentVariability fv;
      fv.setVariability(x, chr.start, chr.end, fwd, rev);
      
      double diffMax(0);
      for(size_t k=0; k<sizeOfvDistOfDistaneOfFrag; ++k) {
	//      std::cout << fv.getAccuOfDistanceOfFragment(k) << "\t" << fvback[k] << std::endl;
	diffMax = std::max(diffMax, fv.getAccuOfDistanceOfFragment(k) - fvback[k]);
      }
      chr.setmp(x, diffMax, mtx);
      mpfv[x].add2genome(fv, mtx);
    }
  }

  return;
}

void scanRepeatRegion(const std::vector<char> &fwd, const std::vector<char> &rev)
{
  int SizeOfFragOverlapDist(10000);
  std::vector<int> FragOverlapDist(SizeOfFragOverlapDist,0);

  int size(fwd.size());
  std::vector<short> array(size, 0);
  int fraglen=1000;
  int last(0);
  for(int i=0; i<size - fraglen; ++i) {
    if(fwd[i] && rev[i+fraglen]) {
      for(int j=0; j<fraglen; ++j) ++array[i+j];
      if(i-last==1) std::cout << last << "-" << i << std::endl;
      last=i;
    }
  }

  //  for(int i=0; i<size; ++i) if(array[i]>1) std::cout << i << "\t" << array[i] << std::endl;
  return;
}
