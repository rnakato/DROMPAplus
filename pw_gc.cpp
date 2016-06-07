#include "pw_gc.h"
#include "pw_makefile.h"


vector<int> readMpbl_binary(const variables_map &values, const SeqStats &chr)
{
  string filename = values["mp"].as<string>() + "/map_" + chr.name + "_binary.txt";
  vector<int> mparray(chr.len, 0);

  isFile(filename);
  int n(0);
  char c;
  ifstream in(filename);
  while (!in.eof()) {
    c = in.get();
    if(c==' ') continue;
    if(c=='1') ++mparray[n];
    ++n;
    if(n >= chr.len-1) break;
  }

  return mparray;
}

/*
vector<int> make_fragarray(const variables_map &values, vector<SeqStats>::iterator chr){
  // mappability
  //  auto ar = readMpbl_binary(values, chr); // base-pair resolution

  // ignore bed regions
    int nbed, beds, bede;
  if(p->bedfilename){
    nbed = p->enrichfile->chr[chr].num;
    for(i=0; i<nbed; i++){
      beds = p->enrichfile->chr[chr].bed[i].s;
      bede = p->enrichfile->chr[chr].bed[i].e;
      for(j=beds; j<=bede; j++) ar[j] = 0;
    }
    }
  return ar;
}*/

  

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

  //  auto ar = make_fragarray(values, p.chr[nchr]);

  /*// make fastaGCarray
  char *fstfile = alloc_str_new(p->genomefile, 50);
  sprintf(fstfile, "%s/%s.fa", p->genomefile, g->chr[chrref].name);

  TYPE_FASTAGCARRAY *fastaGCarray = make_fastaGCarray(fstfile, g->chr[chrref].len, p->flen4gc);
  MYFREE(fstfile);

  // make GCarray for genome
  g->GCdist = makeGCarray_genome(fastaGCarray, ar, g->chr[chrref].len, p->flen4gc);
  // make GCarray for reads
  mapfile->GCdist = makeGCarray_read(p, mapfile, fastaGCarray, ar, chrref, g->chr[chrref].len, p->flen4gc);

  MYFREE(fastaGCarray);
  MYFREE(ar);

  int max=0;
  for(i=0; i<=p->flen4gc; i++){
    g->sum_GCdist       += g->GCdist[i];
    mapfile->sum_GCdist += mapfile->GCdist[i];
    if(max < mapfile->GCdist[i]){
      max = mapfile->GCdist[i];
      mapfile->maxGC = i;
    }
  }
  mapfile->GCweight = (double *)my_calloc(p->flen4gc +1, sizeof(double), "mapfile->GCweight");

  output_GCdist(p, mapfile, g);
  */
  return;
}
