/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "pw_gc.h"
#include "pw_makefile.h"
#include "readdata.h"
#include "macro.h"

#define FRAG_IGNORE 5
#define GCDIST_THRE 1e-5
#define GCDEPTH_THRE 1e-3

vector<short> makeFastaArray(string filename, int length, int flen4gc);
vector<int> makeDistRead(Mapfile &p, vector<short> &fastaGCarray, const vector<char> &ar, SeqStats &chr, int, int, int);
vector<int> makeDistGenome(const vector<short> &FastaArray, const vector<char> &ar, int chrlen, int flen);

void output_GCdist(const variables_map &values, Mapfile &p, const vector<int> &gDist, const vector<int> &rDist, int flen4gc)
{
  int GsumGC = accumulate(gDist.begin(), gDist.end(), 0);
  int RsumGC = accumulate(rDist.begin(), rDist.end(), 0);

  string filename = p.oprefix + ".GCdist.xls";
  ofstream out(filename);
  out << "GC\tgenome prop\treads prop\tdepth\tweight" << endl;
  
  for(int i=0; i<=flen4gc; ++i) {
    double r_genome = gDist[i] / (double)GsumGC;
    double r_read   = rDist[i] / (double)RsumGC;
    double r_depth  = gDist[i] ? rDist[i]/(double)gDist[i]: 0;
    
    out << boost::format("%1%\t%2%\t%3%\t%4%\t") % i % r_genome % r_read % r_depth;
    if(!r_read || r_genome < GCDIST_THRE) p.GCweight[i] = 0;
    else{
      if(!values.count("gcdepthoff") && r_depth < GCDEPTH_THRE) p.GCweight[i] = 1;
      else p.GCweight[i] = r_genome/r_read;
    }
    out << p.GCweight[i] << endl;
  }
  cout << "fragment distribution is output to "<< filename << "." << endl;
  
  return;
}

void make_GCdist(const variables_map &values, Mapfile &p)
{
  cout << "chromosome for GC distribution: chr" << p.lchr->name << endl;
  int chrlen(p.lchr->len);
  int flen(p.getflen(values));
  int flen4gc = min(values["flen4gc"].as<int>(), flen - FRAG_IGNORE*2);
  BPRINT("GC distribution from %1% bp to %2% bp of fragments.\n") % FRAG_IGNORE % (flen4gc+FRAG_IGNORE);

  // mappability
  vector<char> array; 
  if(values.count("mp")) array = readMpbl_binary(values["mp"].as<string>(), ("chr" + p.lchr->name), chrlen);
  else array = readMpbl_binary(p.lchr->len);
  if(values.count("bed")) arraySetBed(array, p.lchr->name, p.vbed);

  // make fastaGCarray
  string fa= values["genome"].as<string>() + "/chr" + p.lchr->name + ".fa";
  auto FastaArray = makeFastaArray(fa, chrlen, flen4gc);

  // make GCarray for genome
  auto gDist = makeDistGenome(FastaArray, array, chrlen, flen4gc);
  // make GCarray for reads
  auto rDist = makeDistRead(p, FastaArray, array, *p.lchr, chrlen, flen, flen4gc);
  p.maxGC = getmaxi(rDist);

  output_GCdist(values, p, gDist, rDist, flen4gc);

  return;
}

vector<int> makeDistRead(Mapfile &p, vector<short> &fastaGCarray, const vector<char> &ar, SeqStats &chr, int chrlen, int flen, int flen4gc)
{
  int gc, posi;

  vector<int> array(flen4gc+1, 0);
  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto x:chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      if(strand==STRAND_PLUS) posi = min(x.F3 + FRAG_IGNORE, chrlen -1);
      else                    posi = max(x.F3 - flen + FRAG_IGNORE, 0);
      if(!ar[posi] || !ar[posi + flen4gc]) continue;
      gc = fastaGCarray[posi];
      if(gc != -1) array[gc]++;
    }
  }

  return array;
}

vector<int> makeDistGenome(const vector<short> &FastaArray, const vector<char> &ar, int chrlen, int flen4gc)
{
  vector<int> array(flen4gc+1, 0);

  int max = chrlen - FRAG_IGNORE - flen4gc;
  for(int i= FRAG_IGNORE + flen4gc; i<max; ++i) {
    if(ar[i]) {
      short gc(FastaArray[i]);
      if(gc != -1) array[gc]++;
    }
  }
  return array;
}


/* return -1 when including Ns */
vector<short> makeFastaArray(string filename, int length, int flen4gc)
{
  int s,e, n(0);
  int state(0);
  char c;
  vector<short> array(length,0);
  ifstream in(filename);
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
	  s = max(0, n - flen4gc);
	  e = n;
	  for(int i=s; i<e; i++){
	    if(array[i] != -1) array[i]++;
	  }
	}else if(c=='A' || c=='T' || c=='a' || c=='t') {
	  /* none */
	}else{  /* N and others */
	  s = max(0, n - flen4gc);
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

void weightReadchr(const variables_map &values, Mapfile &p, int s, int e, boost::mutex &mtx)
{
  for(int i=s; i<=e; ++i) {
    int posi;
    int flen(p.getflen(values));
    int flen4gc = min(values["flen4gc"].as<int>(), flen - FRAG_IGNORE*2);
    cout << p.chr[i].name << ".." << flush;
    string fa = values["genome"].as<string>() + "/chr" + p.chr[i].name + ".fa";
    auto FastaArray = makeFastaArray(fa, p.chr[i].len, flen4gc);
    
    for(int strand=0; strand<STRANDNUM; ++strand) {
      for (auto &x: p.chr[i].seq[strand].vRead) {
	if(x.duplicate) continue;
	if(strand==STRAND_PLUS) posi = min(x.F3 + FRAG_IGNORE, (int)p.chr[i].len -1);
	else                    posi = max(x.F3 - flen + FRAG_IGNORE, 0);
	int gc(FastaArray[posi]);
	if(gc != -1) x.multiplyWeight(p.GCweight[gc]);
	p.chr[i].seq[strand].addReadAfterGC(x.getWeight(), mtx);
      }
    }
  }
  return;
}

void weightRead(const variables_map &values, Mapfile &p)
{
  cout << "add weight to reads..." << flush;

  boost::thread_group agroup;
  boost::mutex mtx;
  for(uint i=0; i<p.vsepchr.size(); i++) {
    agroup.create_thread(bind(weightReadchr, boost::cref(values), boost::ref(p), p.vsepchr[i].s, p.vsepchr[i].e, boost::ref(mtx)));
  }
  agroup.join_all();

  printf("done.\n");

  // Stats
  for(auto &x:p.chr) {
    for(int strand=0; strand<STRANDNUM; ++strand) {
      p.genome.seq[strand].nread_afterGC += x.seq[strand].nread_afterGC;
    }
  }
  return;
}
