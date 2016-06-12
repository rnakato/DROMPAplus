/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */

#include <boost/algorithm/string.hpp>
#include "pw_makefile.h"
#include "readdata.h"

#define PRINTWARNING_W(w) cerr << "Warning: Read scaling weight = " << w << ". Too much scaling up will bring noisy results." << endl;

void addReadToWigArray(const variables_map &, vector<int> &, const Read, long);
vector<int> makeWigarray(const variables_map &, Mapfile &, SeqStats &);
void norm2rpm(const variables_map &, Mapfile &, SeqStats &, vector<int> &);
void outputWig(const variables_map &, Mapfile &, string);
void outputBedGraph(const variables_map &, Mapfile &, string);
void outputBinary(const variables_map &, Mapfile &, string);

void makewig(const variables_map &values, Mapfile &p)
{
  printf("Convert read data to array: \n");
  int oftype = values["of"].as<int>();
  string filename = p.oprefix + "." + IntToString(values["binsize"].as<int>());

  if(oftype==TYPE_COMPRESSWIG || oftype==TYPE_UNCOMPRESSWIG){
    filename += ".wig";
    outputWig(values, p, filename);
    if(oftype==TYPE_COMPRESSWIG) {
      string command = "gzip -f " + filename;
      if(system(command.c_str())) PRINTERR("gzip .wig failed.");
    }
  } else if (oftype==TYPE_BEDGRAPH || oftype==TYPE_BIGWIG) {
    filename += ".bedGraph";
    outputBedGraph(values, p, filename);
    if(oftype==TYPE_BIGWIG) {
      printf("Convert to bigWig...\n");
      string command = "bedGraphToBigWig " + filename + " " + values["gt"].as<string>() + " " + p.oprefix + ".bw";
      if(system(command.c_str())) PRINTERR("conversion failed.");
      remove(filename.c_str()); 
    }
  } else if (oftype==TYPE_BINARY) {
    filename += ".bin";
    outputBinary(values, p, filename);
  }

  printf("done.\n");
  return;
}

void addReadToWigArray(const variables_map &values, vector<int> &wigarray, const Read x, long chrlen)
{
  int s, e;
  s = min(x.F3, x.F5);
  e = max(x.F3, x.F5);
  if(values["rcenter"].as<int>()) {  // consider only center region of fragments
    s = (s + e - values["rcenter"].as<int>())/2;
    e = s + values["rcenter"].as<int>();
  }
  s = max(0, s);
  e = min(e, (int)(chrlen -1));

  int sbin(s/values["binsize"].as<int>());
  int ebin(e/values["binsize"].as<int>());
  for(int j=sbin; j<=ebin; ++j) wigarray[j] += VALUE2WIGARRAY(x.getWeight());
  return;
}

vector<int> makeWigarray(const variables_map &values, Mapfile &p, SeqStats &chr)
{
  vector<int> wigarray(chr.nbin, 0);

  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto x:chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      addReadToWigArray(values, wigarray, x, chr.len);
    }
  }

  // Mappability normalization
  if (values.count("mp")) {
    int binsize = values["binsize"].as<int>();
    int mpthre = values["mpthre"].as<double>()*binsize;
    int mpthre_dist = 0.8*binsize;
    auto mparray = readMpbl(values["mp"].as<string>(), chr.name, values["binsize"].as<int>(), chr.nbin);
    vector<int> mappablearray;
    for(int i=0; i<chr.nbin; ++i) {
      if(mparray[i] > mpthre) wigarray[i] = wigarray[i]*binsize/(double)mparray[i];
      chr.addmpDist(mparray[i]/(double)binsize);
      if(mparray[i] > mpthre_dist) mappablearray.push_back(mparray[i]);
    }
    chr.getWigStats(mappablearray);
  } else {
    chr.getWigStats(wigarray);
  }
  p.genome.addWigStats(chr);
  
  /* Total read normalization */
  if(values["ntype"].as<string>() != "NONE") norm2rpm(values, p, chr, wigarray);

#ifdef DEBUG
  chr.printWigStats();
  chr.printmpDist();
#endif

  return wigarray;
}

template <class T, class S>
double setw(T nm, S dn)
{
  return (dn ? nm/(double)dn: 0);
}

void norm2rpm(const variables_map &values, Mapfile &p, SeqStats &chr, vector<int> &wigarray)
{
  double w(0);
  string ntype(values["ntype"].as<string>());
  
  if(ntype == "GR") {
    double dn(p.genome.bothnread_nonred());
    w = setw(values["nrpm"].as<int>(), dn);
    if(chr.name == p.lastchr) {
      BPRINT("\ngenomic read number = %1%, after=%2%, w=%3$.3f\n") % (long)dn % values["nrpm"].as<int>() % w;
      if(w>2) PRINTWARNING_W(w);
    }
  } else if(ntype == "GD") {
    w = setw(values["ndepth"].as<double>(), p.genome.depth);
    if(chr.name == p.lastchr) {
      BPRINT("\ngenomic depth = %1$.2f, after=%2$.2f, w=%3$.3f\n") % p.genome.depth % values["ndepth"].as<double>() % w;
      if(w>2) PRINTWARNING_W(w);
    }
  } else if(ntype == "CR") {
    double nm = values["nrpm"].as<int>() * (chr.len_mpbl/(double)p.genome.len_mpbl);
    double dn = chr.bothnread_nonred();
    w = setw(nm, dn);
    BPRINT("read number = %1%, after=%2$.1f, w=%3$.3f\n") % (long)dn % nm % w;
    if(w>2) PRINTWARNING_W(w);
  } else if(ntype == "CD") {
    w = setw(values["ndepth"].as<double>(), chr.depth);
    BPRINT("depth = %1$.2f, after=%2$.2f, w=%3$.3f\n") % chr.depth % values["ndepth"].as<double>() % w;
    if(w>2) PRINTWARNING_W(w);
  }

  for(int i=0; i<chr.nbin; ++i) { if(wigarray[i]) wigarray[i] *= w;}

  chr.setWeight(w);
  if(ntype == "GR" || ntype == "GD") p.genome.setWeight(w);
  else { for(int i=0; i<STRANDNUM; i++) p.genome.seq[i].nread_rpm += chr.seq[i].nread_rpm;}
  return;
}

void outputWig(const variables_map &values, Mapfile &p, string filename)
{
  int binsize = values["binsize"].as<int>();
  ofstream out(filename);
  out << boost::format("track type=wiggle_0\tname=\"%1%\"\tdescription=\"Merged tag counts for every %2% bp\"\n")
    % values["output"].as<string>() % binsize;
  
  for(auto &chr: p.chr) {
    vector<int> array = makeWigarray(values, p, chr);
    out << boost::format("variableStep\tchrom=%1%\tspan=%2%\n") % chr.name % binsize;
    for(int i=0; i<chr.nbin; ++i) {
      if(array[i]) out << i*binsize +1 << "\t" << WIGARRAY2VALUE(array[i]) << endl;
    }
  }
  
#ifdef DEBUG
  p.genome.printWigStats();
#endif

  return;
}

void outputBedGraph(const variables_map &values, Mapfile &p, string filename)
{
  int e;
  int binsize = values["binsize"].as<int>();
  
  ofstream out(filename);
  out << boost::format("browser position %1%:%2%-%3%\n") % p.chr[1].name % 0 % (p.chr[1].len/100);
  out << "browser hide all" << endl;
  out << "browser pack refGene encodeRegions" << endl;
  out << "browser full altGraph" << endl;
  out << boost::format("track type=bedGraph name=\"%1%\" description=\"Merged tag counts for every %2% bp\" visibility=full\n")
    % values["output"].as<string>() % binsize;
  out.close();

  string tempfile = filename + ".temp";
  ofstream temp(tempfile);
  for(auto &chr: p.chr) {
    cout << chr.name << ".." << flush;
    vector<int> array = makeWigarray(values, p, chr);
    for(int i=0; i<chr.nbin; ++i) {
      if(i==chr.nbin -1) e = chr.len -1; else e = (i+1)*binsize;
      if(array[i]) temp << chr.name << " " << i*binsize << " " <<e << " " << WIGARRAY2VALUE(array[i]) << endl;
    }
  }
  temp.close();
  
  string command = "sort -k1,1 -k2,2n "+ tempfile +" >> " + filename;
  if(system(command.c_str())) PRINTERR("sorting bedGraph failed.");

  remove(tempfile.c_str());
  
  return;
}

void outputBinary(const variables_map &values, Mapfile &p, string filename)
{
  ofstream out(filename, ios::binary);
  for(auto &chr: p.chr) {
    vector<int> array = makeWigarray(values, p, chr);
    for(int i=0; i<chr.nbin; ++i) out.write((char *)&array[i], sizeof(int));
  }
  return;
}
