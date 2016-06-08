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

vector<int> makeDistRead(Mapfile &p, vector<short> &fastaGCarray, const vector<char> &ar, SeqStats &chr, int chrlen, int flen)
{
  int gc, posi;

  vector<int> array(flen+1, 0);

  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto x:chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      if(strand==STRAND_PLUS) posi = min(x.F3 + FRAG_IGNORE, chrlen -1);
      else                    posi = max(x.F3 - p.dist.eflen + FRAG_IGNORE, 0);
      if(!ar[posi] || !ar[posi + flen]) continue;
      gc = fastaGCarray[posi];
      if(gc != -1) array[gc]++;
    }
  }

  return array;
}

vector<int> makeDistGenome(const vector<short> &FastaArray, const vector<char> &ar, int chrlen, int flen)
{
  vector<int> array(flen+1, 0);

  int max = chrlen - FRAG_IGNORE - flen;
  for(int i= FRAG_IGNORE + flen; i<max; ++i) {
    if(ar[i]) {
      short gc(FastaArray[i]);
      if(gc != -1) array[gc]++;
    }
  }
  return array;
}

void output_GCdist(const variables_map &values, Mapfile &p, const vector<int> &gDist, const vector<int> &rDist, int flen)
{
  int max(0);
  int GsumGC(0),RsumGC(0);
  for(int i=0; i<=flen; ++i){
    GsumGC += gDist[i];
    RsumGC += rDist[i];
    if(max < rDist[i]){
      max = rDist[i];
      p.maxGC = i;
    }
  }
  
  string filename = values["odir"].as<string>() + "/" + values["output"].as<string>() + ".GCdist.xls";
  ofstream out(filename);
  out << "GC\tgenome prop\treads prop\tdepth\tweight" << endl;
  
  for(int i=0; i<=flen; ++i) {
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
  /* use longest chromosome as reference genome */
  long lenmax(0);
  int nchr = 0;
  for(size_t i=0; i<p.chr.size(); ++i) {
    if(lenmax < p.chr[i].len) {
      lenmax = p.chr[i].len;
      nchr = i;
    }
  }
  cout << "chromosome for GC distribution: " << p.chr[nchr].name << endl;
  int chrlen(p.chr[nchr].len);
  int flen(values["flen4gc"].as<int>());
  
  // mappability
  vector<char> array; 
  if(values.count("mp")) array = readMpbl_binary(values["mp"].as<string>(), p.chr[nchr].name, chrlen);
  else array = readMpbl_binary(p.chr[nchr].len);
  if(values.count("bed")) arraySetBed(array, p.chr[nchr].name, p.vbed);
  
  // make fastaGCarray
  string fa= values["genome"].as<string>() + "/" + p.chr[nchr].name + ".fa";
  
  auto FastaArray = makeFastaArray(fa, chrlen, flen);
  
  // make GCarray for genome
  auto gDist = makeDistGenome(FastaArray, array, chrlen, flen);
  // make GCarray for reads
  auto rDist = makeDistRead(p, FastaArray, array, p.chr[nchr], chrlen, flen);
  
  output_GCdist(values, p, gDist, rDist, flen);

  return;
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

void weightRead(const variables_map &values, Mapfile &p)
{
  int gc, posi;
  int flen(values["flen4gc"].as<int>());

  cout << "add weight to reads..." << flush;
  for(auto &x:p.chr) {
    cout << x.name << ".." << flush;
    string fa = values["genome"].as<string>() + "/" + x.name + ".fa";
    auto FastaArray = makeFastaArray(fa, x.len, flen);

    for(int strand=0; strand<STRANDNUM; ++strand) {
      for (auto &read: x.seq[strand].vRead) {
	if(read.duplicate) continue;
	if(strand==STRAND_PLUS) posi = min(read.F3 + FRAG_IGNORE, (int)x.len -1);
	else                    posi = max(read.F3 - p.dist.eflen + FRAG_IGNORE, 0);
	gc = FastaArray[posi];
	if(gc != -1) read.weight *= p.GCweight[gc];
	x.seq[strand].nread_afterGC += read.weight;
      }
    }
  }
  printf("done.\n");

  // Stats
  for(auto &x:p.chr) {
    for(int strand=0; strand<STRANDNUM; ++strand) {
      p.genome.seq[strand].nread_afterGC += x.seq[strand].nread_afterGC;
    }
  }
  return;
}
