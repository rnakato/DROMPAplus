/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "pw_makefile.h"
#include "mytype.h"
#include "readbpstatus.h"

namespace {
  void printwarning(double w)
  {
    std::cerr << "Warning: Read scaling weight = " << w << ". Too much scaling up will bring noisy results." << std::endl;
  }
  
  template <class T, class S>
  double setw(T nm, S dn) { return (dn ? getratio(nm, dn): 0); }
}
  
void addReadToWigArray(const MyOpt::Variables &, WigArray &, const Read, const int64_t);
WigArray makeWigarray(const MyOpt::Variables &, Mapfile &, SeqWigStats &);
void norm2rpm(const MyOpt::Variables &values, Mapfile &p, SeqStats &chr, WigArray &wigarray);
void outputWig(const MyOpt::Variables &, Mapfile &, const std::string &);
void outputBedGraph(const MyOpt::Variables &, Mapfile &, const std::string &);
void outputBinary(const MyOpt::Variables &, Mapfile &, const std::string &);

void makewig(const MyOpt::Variables &values, Mapfile &p)
{
  printf("Convert read data to array: \n");
  WigType oftype(static_cast<WigType>(values["of"].as<int>()));
  std::string filename(p.getbinprefix());

  if (oftype==WigType::COMPRESSWIG || oftype==WigType::UNCOMPRESSWIG) {
    filename += ".wig";
    outputWig(values, p, filename);
    if(oftype==WigType::COMPRESSWIG) {
      std::string command = "gzip -f " + filename;
      if(system(command.c_str())) PRINTERR("gzip .wig failed.");
    }
  } else if (oftype==WigType::BEDGRAPH || oftype==WigType::BIGWIG) {
    filename += ".bedGraph";
    outputBedGraph(values, p, filename);
    if(oftype==WigType::BIGWIG) {
      printf("Convert to bigWig...\n");
      std::string command = "bedGraphToBigWig " + filename + " " + values["gt"].as<std::string>() + " " + p.getprefix() + ".bw";
      if(system(command.c_str())) PRINTERR("conversion failed.");
      remove(filename.c_str()); 
    }
  } else if (oftype==WigType::BINARY) {
    filename += ".bin";
    outputBinary(values, p, filename);
  }

  printf("done.\n");
  return;
}

void addReadToWigArray(const MyOpt::Variables &values, WigArray &wigarray, const Read x, const int64_t chrlen)
{
  int s, e;
  s = std::min(x.F3, x.F5);
  e = std::max(x.F3, x.F5);

  int rcenter(values["rcenter"].as<int>());
  if(rcenter) {  // consider only center region of fragments
    s = (s + e - rcenter)/2;
    e = s + rcenter;
  }
  s = std::max(0, s);
  e = std::min(e, (int)(chrlen -1));

  int sbin(s/values["binsize"].as<int>());
  int ebin(e/values["binsize"].as<int>());
  for(int j=sbin; j<=ebin; ++j) wigarray.addval(j, x.getWeight());
  return;
}

void peakcall(Mapfile &mapfile, const SeqWigStats &chr, const WigArray &wigarray)
{
  int size = wigarray.size();
  int ext(0);
  double pthre(3);

  for(int i=0; i<size; ++i) {
    double val(wigarray.getval(i));
    double p(getlogpZINB(val, chr.ws.nb_p, chr.ws.nb_n));

    if(!ext) {
      if(p > pthre) {
	mapfile.addPeak(Peak(i, i, chr.getname(), val, p));
	ext=1;
      }
    } else {
      if(p > pthre) {
	mapfile.renewPeak(i, val, p);
      } else {
	ext=0;
      }
    }
  }
  return;
}

WigArray makeWigarray(const MyOpt::Variables &values, Mapfile &p, SeqWigStats &chr)
{
  std::cout << chr.getname() << ".." << std::flush;
  WigArray wigarray(chr.getnbin(), 0);

  for (auto strand: {Strand::FWD, Strand::REV}) {
    for (auto &x: chr.getvReadref(strand)) {
      if(x.duplicate) continue;
      addReadToWigArray(values, wigarray, x, chr.getlen());
    }
  }

  if (values.count("mp")) {
    int binsize = values["binsize"].as<int>();
    int mpthre = values["mpthre"].as<double>()*binsize;
    auto mparray = readMpbl(values["mp"].as<std::string>(), ("chr" + chr.getname()), values["binsize"].as<int>(), chr.getnbin());
    for(int i=0; i<chr.getnbin(); ++i) {
      chr.ws.addmpDist(getratio(mparray[i], binsize));
      if(mparray[i] > mpthre) wigarray.multipleval(i, getratio(binsize, mparray[i]));
    }
  }
  chr.ws.getWigStats(wigarray);
  p.wsGenome.addWigDist(chr.ws);

  peakcall(p, chr, wigarray);
  /* Total read normalization */
  if(values["ntype"].as<std::string>() != "NONE") norm2rpm(values, p, chr, wigarray);

#ifdef DEBUG
  if (values.count("mp")) chr.ws.printmpDist();
#endif

  return wigarray;
}

void norm2rpm(const MyOpt::Variables &values, Mapfile &p, SeqStats &chr, WigArray &wigarray)
{
  static int on(0);
  double w(0);
  std::string ntype(values["ntype"].as<std::string>());
  
  if(ntype == "GR") {
    double dn(p.genome.getnread_nonred(Strand::BOTH));
    w = setw(values["nrpm"].as<int>(), dn);
    if(!on) {
      std::cout << boost::format("\ngenomic read number = %1%, after=%2%, w=%3$.3f\n") % (int64_t)dn % values["nrpm"].as<int>() % w;
      if(w>2) printwarning(w);
      on=1;
    }
  } else if(ntype == "GD") {
    w = setw(values["ndepth"].as<double>(), p.genome.getdepth());
    if(!on) {
      std::cout << boost::format("\ngenomic depth = %1$.2f, after=%2$.2f, w=%3$.3f\n") % p.genome.getdepth() % values["ndepth"].as<double>() % w;
      if(w>2) printwarning(w);
      on=1;
    }
  } else if(ntype == "CR") {
    double nm = values["nrpm"].as<int>() * getratio(chr.getlenmpbl(), p.genome.getlenmpbl());
    double dn = chr.getnread_nonred(Strand::BOTH);
    w = setw(nm, dn);
    std::cout << boost::format("read number = %1%, after=%2$.1f, w=%3$.3f\n") % static_cast<int64_t>(dn) % nm % w;
    if(w>2) printwarning(w);
  } else if(ntype == "CD") {
    w = setw(values["ndepth"].as<double>(), chr.getdepth());
    std::cout << boost::format("depth = %1$.2f, after=%2$.2f, w=%3$.3f\n") % chr.getdepth() % values["ndepth"].as<double>() % w;
    if(w>2) printwarning(w);
  }

  for(int i=0; i<chr.getnbin(); ++i) { wigarray.multipleval(i, w);}

  chr.setsizefactor(w);
  if(ntype == "GR" || ntype == "GD") p.genome.setsizefactor(w);
  return;
}

void outputWig(const MyOpt::Variables &values, Mapfile &p, const std::string &filename)
{
  int binsize = values["binsize"].as<int>();
  std::ofstream out(filename);
  out << boost::format("track type=wiggle_0\tname=\"%1%\"\tdescription=\"Merged tag counts for every %2% bp\"\n")
    % values["output"].as<std::string>() % binsize;
  
  for(auto &chr: p.genome.chr) {
    WigArray array = makeWigarray(values, p, chr);
    out << boost::format("variableStep\tchrom=%1%\tspan=%2%\n") % chr.getname() % binsize;
    array.outputAsWig(out, binsize);
  }
  
  return;
}

void outputBedGraph(const MyOpt::Variables &values, Mapfile &p, const std::string &filename)
{
  int32_t binsize = values["binsize"].as<int>();
  
  std::ofstream out(filename);
  out << boost::format("browser position %1%:%2%-%3%\n") % p.genome.chr[1].getname() % 0 % (p.genome.chr[1].getlen()/100);
  out << "browser hide all" << std::endl;
  out << "browser pack refGene encodeRegions" << std::endl;
  out << "browser full altGraph" << std::endl;
  out << boost::format("track type=bedGraph name=\"%1%\" description=\"Merged tag counts for every %2% bp\" visibility=full\n")
    % values["output"].as<std::string>() % binsize;
  out.close();

  std::string tempfile = filename + ".temp";
  std::ofstream out2(tempfile);
  for(auto &chr: p.genome.chr) {
    WigArray array = makeWigarray(values, p, chr);
    array.outputAsBedGraph(out2, binsize, chr.getname(), chr.getlen() -1);
    
    /*    for(int i=0; i<chr.getnbin(); ++i) {
      if(i==chr.getnbin() -1) e = chr.getlen() -1; else e = (i+1)*binsize;
      if(array[i]) out2 << chr.getname() << " " << i*binsize << " " << e << " " << wigarray2value(array[i]) << std::endl;
      }*/
  }
  out2.close();
  
  std::string command = "sort -k1,1 -k2,2n "+ tempfile +" >> " + filename;
  if(system(command.c_str())) PRINTERR("sorting bedGraph failed.");

  remove(tempfile.c_str());
  
  return;
}

void outputBinary(const MyOpt::Variables &values, Mapfile &p, const std::string &filename)
{
  std::ofstream out(filename, std::ios::binary);
  for(auto &chr: p.genome.chr) {
    WigArray array = makeWigarray(values, p, chr);
    array.outputAsBinary(out);
  }
  return;
}
