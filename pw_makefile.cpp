/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */

#include "pw_makefile.h"
#include "macro.h"
using namespace boost::program_options;

vector<int> makeWigarray(const variables_map &, SeqStats &);
void outputWig(const variables_map &, Mapfile &, string);

void makewig(const variables_map &values, Mapfile &p)
{
  printf("Convert read data to array: \n");

  /*  p.wstats.thre = define_wstats_thre(p, mapfile, g);
  p.wstats.n_darray = max(NUM_DARRAY, p.wstats.thre);
  p.wstats.genome->darray_all = (int *)my_calloc(p.wstats.n_darray +1, sizeof(int), "wstats.genome->darray_all");
  p.wstats.genome->darray_bg  = (int *)my_calloc(p.wstats.n_darray +1, sizeof(int), "wstats.genome->darray_bg");
  for(chr=1; chr<g->chrnum; chr++){
    p.wstats.chr[chr].darray_all = (int *)my_calloc(p.wstats.n_darray +1, sizeof(int), "wstats.chr[chr]->darray_all");
    p.wstats.chr[chr].darray_bg  = (int *)my_calloc(p.wstats.n_darray +1, sizeof(int), "wstats.chr[chr]->darray_bg");
    }*/

  int oftype = values["of"].as<int>();
  string filename = values["odir"].as<string>() + "/" + values["output"].as<string>();

  if(oftype==TYPE_COMPRESSWIG || oftype==TYPE_UNCOMPRESSWIG){
    filename += ".wig";
    outputWig(values, p, filename);
    if(oftype==TYPE_COMPRESSWIG) {
      string command = "gzip -f " + filename;
      if(system(command.c_str())) PRINTERR("gzip .wig failed.");
    }
  } else if (oftype==TYPE_BEDGRAPH || oftype==TYPE_BIGWIG) {
    filename += ".bedGraph";
    ofstream out(filename);
    /*    for(auto chr: p.chr) {
      vector<int> array = makeWigarray(values, chr);
      outputBedGraph(values, out, array, chr);
    }
    out.close();
    if(oftype==TYPE_BIGWIG) convert_bedGraph_to_bigWig(outputfilename, output_prefix, gtfile);*/
  } else if (oftype==TYPE_BINARY) {
  }
  
  printf("done.\n");
  return;
}

void addReadToWigArray(const variables_map &values, vector<int> &wigarray, const Read x, long chrlen)
{
  int s, e;
  s = min(x.F3, x.F5);
  e = max(x.F3, x.F5);
  if(values["rcenter"].as<int>()) {  /* consider only center region of fragments */
    s = (s + e - values["rcenter"].as<int>())/2;
    e = s + values["rcenter"].as<int>();
  }
  s = max(0, s);
  e = min(e, (int)(chrlen -1));

  int sbin(s/values["binsize"].as<int>());
  int ebin(e/values["binsize"].as<int>());
  for(int j=sbin; j<=ebin; ++j) wigarray[j] += 1; //VALUE2WIGARRAY(1);
  return;
}

vector<int> makeWigarray(const variables_map &values, SeqStats &chr)
{
  vector<int> wigarray(chr.nbin, 0);

  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto x:chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      addReadToWigArray(values, wigarray, x, chr.len);
    }
  }
  return wigarray;
}

void outputWig(const variables_map &values, Mapfile &p, string filename)
{
  int binsize = values["binsize"].as<int>();
  ofstream out(filename);
  out << boost::format("track type=wiggle_0\tname=\"%1%\"\tdescription=\"Merged tag counts for every %2% bp\"\n")
    % values["output"].as<string>() % binsize;
  
  for(auto chr: p.chr) {
    vector<int> array = makeWigarray(values, chr);
    out << boost::format("variableStep\tchrom=%1%\tspan=%2%\n") % chr.name % binsize;
    for(int i=0; i<chr.nbin; ++i) {
      if(array[i]) out << boost::format("%1%\t%2%\n") % (i*binsize +1) % (double)array[i];
    }
  }
  
  return;
}
