/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "pw_gc.h"

namespace {
  const int lenIgnoreOfFragment(5);
  const double threGcDist(1e-5);
  const double threGcDepth(1e-3);
}

std::vector<int> makeDistGenome(const std::vector<short> &, const std::vector<char> &, const int, const int);
std::vector<int> makeDistRead(const std::vector<short> &, const std::vector<char> &, const SeqStats &, const int, const int, const int);
std::vector<short> makeFastaArray(const std::string &, const int, const int);

class GCdist {
  int flen;
  int flen4gc;
  std::vector<int> DistGenome;
  std::vector<int> DistRead;
  std::vector<double> GCweight;

  double getPropGenome(const int i) {
    double GsumGC = accumulate(DistGenome.begin(), DistGenome.end(), 0);
    return GsumGC ? DistGenome[i] / GsumGC: 0;
  }
  double getPropRead(const int i) {
    double RsumGC = accumulate(DistRead.begin(), DistRead.end(), 0);
    return RsumGC ? DistRead[i] / RsumGC: 0;
  }
  double getPropDepth(const int i) {
    return DistGenome[i] ? DistRead[i]/static_cast<double>(DistGenome[i]): 0;
  }

  void makeGCweightDist(const MyOpt::Variables &);

public:
  GCdist(const MyOpt::Variables &, const Mapfile &);
  int getmaxGC() const { return getmaxi(DistRead); }
  double getGCweight(const int i) const { return GCweight[i]; }
  void outputGCweightDist(const std::string &filename) {
    std::ofstream out(filename);
    out << "GC\tgenome prop\treads prop\tdepth\tweight" << std::endl;
    for(int i=0; i<=flen4gc; ++i) {
      out << boost::format("%1%\t%2%\t%3%\t%4%\t%5%\n")
	% i % getPropGenome(i) % getPropRead(i) % getPropDepth(i) % GCweight[i];
    }
    std::cout << "fragment distribution is output to "<< filename << "." << std::endl;
  }
};

GCdist::GCdist(const MyOpt::Variables &values, const Mapfile &p)
  : flen(p.getflen(values))
{
  flen4gc = std::min(values["flen4gc"].as<int>(), flen - lenIgnoreOfFragment*2);
  BPRINT("GC distribution from %1% bp to %2% bp of fragments.\n") % lenIgnoreOfFragment % (flen4gc + lenIgnoreOfFragment);
  
  int chrlen(p.lchr->getlen());
  std::string chrname(p.lchr->name);
  std::vector<char> mparray; 
  if(values.count("mp")) mparray = readMpbl_binary(values["mp"].as<std::string>(), ("chr" + chrname), chrlen);
  else mparray = readMpbl_binary(chrlen);
  if(values.count("bed")) arraySetBed(mparray, chrname, p.genome.getvbed());
  
  std::string fastaname= values["genome"].as<std::string>() + "/chr" + chrname + ".fa";
  auto FastaArray = makeFastaArray(fastaname, chrlen, flen4gc);
  
  DistGenome = makeDistGenome(FastaArray, mparray, chrlen, flen4gc);
  DistRead = makeDistRead(FastaArray, mparray, *p.lchr, chrlen, flen, flen4gc);

  makeGCweightDist(values);
}

void GCdist::makeGCweightDist(const MyOpt::Variables &values)
{
  for(int i=0; i<=flen4gc; ++i) {
    double r_genome = getPropGenome(i);
    double r_read   = getPropRead(i);
    double r_depth  = getPropDepth(i);
    
    double w(0);
    if(!r_read || r_genome < threGcDist) w = 0;
    else{
      if(!values.count("gcdepthoff") && r_depth < threGcDepth) w = 1;
      else w = r_genome/r_read;
    }
    GCweight.push_back(w);
  }
}

void weightReadchr(const MyOpt::Variables &values, Mapfile &p, GCdist &dist, int s, int e, boost::mutex &mtx)
{
  for(int i=s; i<=e; ++i) {
    int posi;
    int flen(p.getflen(values));
    int flen4gc = std::min(values["flen4gc"].as<int>(), flen - lenIgnoreOfFragment*2);
    std::cout << p.genome.chr[i].name << ".." << std::flush;
    std::string fa = values["genome"].as<std::string>() + "/chr" + p.genome.chr[i].name + ".fa";
    auto FastaArray = makeFastaArray(fa, p.genome.chr[i].getlen(), flen4gc);
    
    for(int strand=0; strand<STRANDNUM; ++strand) {
      for (auto &x: p.genome.chr[i].seq[strand].vRead) {
	if(x.duplicate) continue;
	if(strand==STRAND_PLUS) posi = std::min(x.F3 + lenIgnoreOfFragment, (int)p.genome.chr[i].getlen() -1);
	else                    posi = std::max(x.F3 - flen + lenIgnoreOfFragment, 0);
	int gc(FastaArray[posi]);
	if(gc != -1) x.multiplyWeight(dist.getGCweight(gc));
	p.genome.chr[i].seq[strand].addReadAfterGC(x.getWeight(), mtx);
      }
    }
  }
  return;
}

void weightRead(const MyOpt::Variables &values, Mapfile &p, GCdist &dist)
{
  std::cout << "add weight to reads..." << std::flush;

  boost::thread_group agroup;
  boost::mutex mtx;
  for(uint i=0; i<p.genome.vsepchr.size(); i++) {
    agroup.create_thread(bind(weightReadchr, boost::cref(values), boost::ref(p), boost::ref(dist), p.genome.vsepchr[i].s, p.genome.vsepchr[i].e, boost::ref(mtx)));
  }
  agroup.join_all();

  printf("done.\n");

  // Stats
  for(auto &x:p.genome.chr) {
    for(int strand=0; strand<STRANDNUM; ++strand) {
      p.genome.seq[strand].nread_afterGC += x.seq[strand].nread_afterGC;
    }
  }
  return;
}

void normalizeByGCcontents(const MyOpt::Variables &values, Mapfile &p)
{
  std::cout << "chromosome for GC distribution: chr" << p.lchr->name << std::endl;
  GCdist dist(values, p);

  p.setmaxGC(dist.getmaxGC());

  std::string filename = p.getprefix() + ".GCdist.xls";
  dist.outputGCweightDist(filename);

  weightRead(values, p, dist);

  return;
}

std::vector<int> makeDistGenome(const std::vector<short> &FastaArray, const std::vector<char> &mparray, const int chrlen, const int flen4gc)
{
  std::vector<int> array(flen4gc+1, 0);

  int end = chrlen - lenIgnoreOfFragment - flen4gc;
  for(int i= lenIgnoreOfFragment + flen4gc; i<end; ++i) {
    if(mparray[i]) {
      short gc(FastaArray[i]);
      if(gc != -1) array[gc]++;
    }
  }
  return array;
}

std::vector<int> makeDistRead(const std::vector<short> &fastaGCarray, const std::vector<char> &mparray, const SeqStats &chr, const int chrlen, const int flen, const int flen4gc)
{
  int posi;
  std::vector<int> array(flen4gc+1, 0);
  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto x:chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      if(strand==STRAND_PLUS) posi = std::min(x.F3 + lenIgnoreOfFragment, chrlen -1);
      else                    posi = std::max(x.F3 - flen + lenIgnoreOfFragment, 0);
      if(!mparray[posi] || !mparray[posi + flen4gc]) continue;
      int gc = fastaGCarray[posi];
      if(gc != -1) array[gc]++;
    }
  }
  return array;
}

/* return -1 when including Ns */
std::vector<short> makeFastaArray(const std::string &filename, const int length, const int flen4gc)
{
  int s,e, n(0);
  int state(0);
  char c;
  std::vector<short> array(length,0);
  std::ifstream in(filename);
  if(!in) PRINTERR("Could not open " << filename << ".");
  
  while (!in.eof()) {
    c = in.get();
    switch(state) {
    case 0:
      if(c=='>') state=1;
      break;
    case 1:    /*header*/
      if(c=='\n') state=2;
      break;
    case 2:    /*body*/
      if(c=='>') goto final;
      else if(isalpha((int)c)) {
	if(c=='G' || c=='C' || c=='g' || c=='c') {
	  s = std::max(0, n - flen4gc);
	  e = n;
	  for(int i=s; i<e; i++){
	    if(array[i] != -1) array[i]++;
	  }
	}else if(c=='A' || c=='T' || c=='a' || c=='t') {
	  /* none */
	}else{  /* N and others */
	  s = std::max(0, n - flen4gc);
	  e = n;
	  for(int i=s; i<e; i++) array[i] = -1;
	}

	n++;
	if(n > length) PRINTERR("ERROR: length " << length << " < " << n);
      }
      break;
    }
  }

 final:
  return array;
}
