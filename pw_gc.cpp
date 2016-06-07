/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "pw_gc.h"
#include "pw_makefile.h"
#include "readdata.h"

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

  // mappability
  vector<char> array; 
  if(values.count("mp")) array = readMpbl_binary(values["mp"].as<string>(), p.chr[nchr].name, p.chr[nchr].len);
  else array = readMpbl_binary(p.chr[nchr].len);
  if(values.count("bed")) arraySetBed(array, p.chr[nchr].name, p.vbed);
  
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
