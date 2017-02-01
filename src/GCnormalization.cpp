/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "GCnormalization.hpp"
#include "WigStats.hpp"
#include "SSP/src/Mapfile.hpp"
#include "SSP/common/util.hpp"

namespace {
  const int32_t lenIgnoreOfFragment=5;
  const double threGcDist(1e-5);
  const double threGcDepth(1e-3);

  std::vector<int> makeDistGenome(const std::vector<short> &, const std::vector<BpStatus> &, const int, const int);
  std::vector<int> makeDistRead(const std::vector<short> &, const std::vector<BpStatus> &, const SeqStats &, const int, const int, const int);
  std::vector<short> makeFastaArray(const std::string &, const int, const int);

  std::vector<int> makeDistGenome(const std::vector<short> &FastaArray, const std::vector<BpStatus> &mparray, const int chrlen, const int flen4gc)
  {
    std::vector<int> array(flen4gc+1, 0);

    int end = chrlen - lenIgnoreOfFragment - flen4gc;
    for(int i= lenIgnoreOfFragment + flen4gc; i<end; ++i) {
      if(mparray[i] != BpStatus::UNMAPPABLE) {
	short gc(FastaArray[i]);
	if(gc != -1) array[gc]++;
      }
    }
    return array;
  }

  std::vector<int> makeDistRead(const std::vector<short> &fastaGCarray, const std::vector<BpStatus> &mparray, const SeqStats &chr, const int chrlen, const int flen, const int flen4gc)
  {
    int posi;
    std::vector<int> array(flen4gc+1, 0);
    for (auto strand: {Strand::FWD, Strand::REV}) {
      for (auto &x: chr.getvReadref(strand)) {
	if(x.duplicate) continue;
	if(strand==Strand::FWD) posi = std::min(x.F3 + lenIgnoreOfFragment, chrlen -1);
	else                    posi = std::max(x.F3 - flen + lenIgnoreOfFragment, 0);
	if(mparray[posi]==BpStatus::UNMAPPABLE || mparray[posi + flen4gc]==BpStatus::UNMAPPABLE) continue;
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
}


  GCdist::GCdist(const int32_t l, GCnorm &gc): flen(l)
  {
    flen4gc = std::min(gc.getflen4gc(), flen - lenIgnoreOfFragment*2);
    std::cout << boost::format("GC distribution from %1% bp to %2% bp of fragments.\n") % lenIgnoreOfFragment % (flen4gc + lenIgnoreOfFragment);
  }
  
  void GCdist::calcGCdist(const SeqStats &chr, const GCnorm &gc, const std::string &mpdir, const int32_t isBedOn, const std::vector<bed> &vbed)
  {
    std::vector<BpStatus> mparray; 
    if(mpdir != "") mparray = readMpbl_binary(mpdir, ("chr" + chr.getname()), chr.getlen());
    else            mparray = readMpbl_binary(chr.getlen());
    if(isBedOn) OverrideBedToArray(mparray, chr.getname(), vbed);
  
    std::string fastaname= gc.getGCdir() + "/chr" + chr.getname() + ".fa";
    auto FastaArray = makeFastaArray(fastaname, chr.getlen(), flen4gc);
  
    DistGenome = makeDistGenome(FastaArray, mparray, chr.getlen(), flen4gc);
    DistRead = makeDistRead(FastaArray, mparray, chr, chr.getlen(), flen, flen4gc);

    makeGCweightDist(gc.isGcDepthOff());
  }
  
  
  void GCdist::makeGCweightDist(const int32_t gcdepthoff)
  {
    for(int i=0; i<=flen4gc; ++i) {
      double r_genome = getPropGenome(i);
      double r_read   = getPropRead(i);
      double r_depth  = getPropDepth(i);
    
      double w(0);
      if(!r_read || r_genome < threGcDist) w = 0;
      else{
	if(!gcdepthoff && r_depth < threGcDepth) w = 1;
	else w = r_genome/r_read;
      }
      GCweight.push_back(w);
    }
  }

  void GCdist::outputGCweightDist(const std::string &filename) {
    std::ofstream out(filename);
    out << "GC\tgenome prop\treads prop\tdepth\tweight" << std::endl;
    for(int32_t i=0; i<=flen4gc; ++i) {
      out << boost::format("%1%\t%2%\t%3%\t%4%\t%5%\n")
	% i % getPropGenome(i) % getPropRead(i) % getPropDepth(i) % GCweight[i];
    }
    std::cout << "fragment distribution is output to "<< filename << "." << std::endl;
  }
  
  void weightReadchr(SeqStatsGenome &genome, GCdist &dist, const std::string &GCdir, int s, int e, boost::mutex &mtx)
  {
    int flen(dist.getflen());
    for(int i=s; i<=e; ++i) {
      int posi;
      std::cout << genome.chr[i].getname() << ".." << std::flush;
      std::string fa = GCdir + "/chr" + genome.chr[i].getname() + ".fa";
      auto FastaArray = makeFastaArray(fa, genome.chr[i].getlen(), dist.getflen4gc());
    
      for (auto strand: {Strand::FWD, Strand::REV}) {
	for (auto &x: genome.chr[i].getvReadref_notconst(strand)) {
	  if(x.duplicate) continue;
	  if(strand==Strand::FWD) posi = std::min(x.F3 + lenIgnoreOfFragment, (int)genome.chr[i].getlen() -1);
	  else                    posi = std::max(x.F3 - flen + lenIgnoreOfFragment, 0);
	  int gc(FastaArray[posi]);
	  if(gc != -1) x.multiplyWeight(dist.getGCweight(gc));

	  genome.chr[i].addReadAfterGC(strand, x.getWeight(), mtx);
	}
      }
    }
    return;
  }


void weightRead(SeqStatsGenome &genome, GCdist &dist, const std::string &GCdir)
{
  std::cout << "add weight to reads..." << std::flush;

  boost::thread_group agroup;
  boost::mutex mtx;
  for(uint i=0; i<genome.vsepchr.size(); i++) {
    agroup.create_thread(bind(weightReadchr, boost::ref(genome), boost::ref(dist), boost::cref(GCdir), genome.vsepchr[i].s, genome.vsepchr[i].e, boost::ref(mtx)));
  }
  agroup.join_all();

  printf("done.\n");
  return;
}

