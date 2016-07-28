/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "pw_gc.h"

namespace {
  const int lenIgnoreOfFragment(5);
  const double threGcDist(1e-5);
  const double threGcDepth(1e-3);
}

std::vector<short> makeFastaArray(std::string filename, int length, int flen4gc);
std::vector<int> makeDistRead(Mapfile &p, std::vector<short> &fastaGCarray, const std::vector<char> &ar, SeqStats &chr, int, int, int);
std::vector<int> makeDistGenome(const std::vector<short> &FastaArray, const std::vector<char> &ar, int chrlen, int flen);

void output_GCdist(const MyOpt::Variables &values, Mapfile &p, const std::vector<int> &gDist, const std::vector<int> &rDist, int flen4gc)
{
  int GsumGC = accumulate(gDist.begin(), gDist.end(), 0);
  int RsumGC = accumulate(rDist.begin(), rDist.end(), 0);

  std::string filename = p.getprefix() + ".GCdist.xls";
  std::ofstream out(filename);
  out << "GC\tgenome prop\treads prop\tdepth\tweight" << std::endl;
  
  for(int i=0; i<=flen4gc; ++i) {
    double r_genome = gDist[i] / static_cast<double>(GsumGC);
    double r_read   = rDist[i] / static_cast<double>(RsumGC);
    double r_depth  = gDist[i] ? rDist[i]/static_cast<double>(gDist[i]): 0;
    
    out << boost::format("%1%\t%2%\t%3%\t%4%\t") % i % r_genome % r_read % r_depth;
    if(!r_read || r_genome < threGcDist) p.GCweight[i] = 0;
    else{
      if(!values.count("gcdepthoff") && r_depth < threGcDepth) p.GCweight[i] = 1;
      else p.GCweight[i] = r_genome/r_read;
    }
    out << p.GCweight[i] << std::endl;
  }
  std::cout << "fragment distribution is output to "<< filename << "." << std::endl;
  
  return;
}

void make_GCdist(const MyOpt::Variables &values, Mapfile &p)
{
  std::cout << "chromosome for GC distribution: chr" << p.lchr->name << std::endl;
  int chrlen(p.lchr->getlen());
  int flen(p.getflen(values));
  int flen4gc = std::min(values["flen4gc"].as<int>(), flen - lenIgnoreOfFragment*2);
  BPRINT("GC distribution from %1% bp to %2% bp of fragments.\n") % lenIgnoreOfFragment % (flen4gc+lenIgnoreOfFragment);

  // mappability
  std::vector<char> array; 
  if(values.count("mp")) array = readMpbl_binary(values["mp"].as<std::string>(), ("chr" + p.lchr->name), chrlen);
  else array = readMpbl_binary(p.lchr->getlen());
  if(values.count("bed")) arraySetBed(array, p.lchr->name, p.genome.getvbed());

  // make fastaGCarray
  std::string fa= values["genome"].as<std::string>() + "/chr" + p.lchr->name + ".fa";
  auto FastaArray = makeFastaArray(fa, chrlen, flen4gc);

  // make GCarray for genome
  auto gDist = makeDistGenome(FastaArray, array, chrlen, flen4gc);
  // make GCarray for reads
  auto rDist = makeDistRead(p, FastaArray, array, *p.lchr, chrlen, flen, flen4gc);
  p.maxGC = getmaxi(rDist);

  output_GCdist(values, p, gDist, rDist, flen4gc);

  return;
}

std::vector<int> makeDistRead(Mapfile &p, std::vector<short> &fastaGCarray, const std::vector<char> &ar, SeqStats &chr, int chrlen, int flen, int flen4gc)
{
  int gc, posi;

  std::vector<int> array(flen4gc+1, 0);
  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto x:chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      if(strand==STRAND_PLUS) posi = std::min(x.F3 + lenIgnoreOfFragment, chrlen -1);
      else                    posi = std::max(x.F3 - flen + lenIgnoreOfFragment, 0);
      if(!ar[posi] || !ar[posi + flen4gc]) continue;
      gc = fastaGCarray[posi];
      if(gc != -1) array[gc]++;
    }
  }

  return array;
}

std::vector<int> makeDistGenome(const std::vector<short> &FastaArray, const std::vector<char> &ar, int chrlen, int flen4gc)
{
  std::vector<int> array(flen4gc+1, 0);

  int max = chrlen - lenIgnoreOfFragment - flen4gc;
  for(int i= lenIgnoreOfFragment + flen4gc; i<max; ++i) {
    if(ar[i]) {
      short gc(FastaArray[i]);
      if(gc != -1) array[gc]++;
    }
  }
  return array;
}


/* return -1 when including Ns */
std::vector<short> makeFastaArray(std::string filename, int length, int flen4gc)
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

void weightReadchr(const MyOpt::Variables &values, Mapfile &p, int s, int e, boost::mutex &mtx)
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
	if(gc != -1) x.multiplyWeight(p.GCweight[gc]);
	p.genome.chr[i].seq[strand].addReadAfterGC(x.getWeight(), mtx);
      }
    }
  }
  return;
}

void weightRead(const MyOpt::Variables &values, Mapfile &p)
{
  std::cout << "add weight to reads..." << std::flush;

  boost::thread_group agroup;
  boost::mutex mtx;
  for(uint i=0; i<p.genome.vsepchr.size(); i++) {
    agroup.create_thread(bind(weightReadchr, boost::cref(values), boost::ref(p), p.genome.vsepchr[i].s, p.genome.vsepchr[i].e, boost::ref(mtx)));
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
