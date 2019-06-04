/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "pw_makefile.hpp"
#include "pw_gv.hpp"
#include "WigStats.hpp"
#include "SSP/src/ReadBpStatus.hpp"
#include "SSP/src/SeqStats.hpp"
#include "SSP/common/inline.hpp"

namespace {
  void printwarning(double w)
  {
    std::cerr << "Warning: Read scaling weight = "
              << w
              << ". Too much scaling up will bring noisy results."
              << std::endl;
  }
}

WigArray makeWigarray(Mapfile &, int32_t);
void norm2rpm(Mapfile &p, SeqStats &chr, WigArray &wigarray);
void outputWig(Mapfile &, const std::string &);
void outputBedGraph(Mapfile &, const std::string &);
//void outputBinary(Mapfile &, const std::string &);

void makewig(Mapfile &p)
{
  printf("Convert read data to array: \n");
  WigType oftype(p.wsGenome.getWigType());
  std::string filename(p.getbinprefix());

  if (oftype==WigType::COMPRESSWIG || oftype==WigType::UNCOMPRESSWIG) {
    filename += ".wig";
    outputWig(p, filename);
    if(oftype==WigType::COMPRESSWIG) {
      std::string command = "gzip -f " + filename;
      if(system(command.c_str())) PRINTERR("gzip .wig failed.");
    }
  } else if (oftype==WigType::BEDGRAPH) {
    filename += ".bedGraph";
    outputBedGraph(p, filename);
  } else if (oftype==WigType::BIGWIG) {
    int32_t fd(0);
    char tmpfile[] = "/tmp/parse2wig+_bedGraph_XXXXXX";
    if ((fd = mkstemp(tmpfile)) < 0){
      perror("mkstemp");
      //      fd = EXIT_FAILURE;
    }
    //    std::cout << "Temporal bedGraph file: " << std::string(tmpfile) << std::endl;
    outputBedGraph(p, std::string(tmpfile));
    printf("Convert to bigWig...\n");
    std::string command = "bedGraphToBigWig " + std::string(tmpfile) + " " + p.genome.getGenomeTable() + " " + p.getbinprefix() + ".bw";
    if(system(command.c_str())) {
      unlink(tmpfile);
      PRINTERR("conversion failed.");
    }
    unlink(tmpfile);
  } /*else if (oftype==WigType::BINARY) {
    filename += ".bin";
    outputBinary(p, filename);
    }*/

  printf("done.\n");
  return;
}

void addReadToWigArray(const WigStatsGenome &p, WigArray &wigarray, const Read x, const int64_t chrlen)
{
  int32_t s, e;
  s = std::min(x.F3, x.F5);
  e = std::max(x.F3, x.F5);

  int32_t rcenter(p.getrcenter());
  if(rcenter) {  // consider only center region of fragments
    s = (s + e - rcenter)/2;
    e = s + rcenter;
  }
  s = std::max(0, s);
  e = std::min(e, (int)(chrlen -1));

  int sbin(s/p.getbinsize());
  int ebin(e/p.getbinsize());
  for(int j=sbin; j<=ebin; ++j) wigarray.addval(j, x.getWeight());
  return;
}

WigArray makeWigarray(Mapfile &p, const int32_t id)
{
  std::cout << "chr" << p.genome.chr[id].getname() << ".." << std::flush;
  WigArray wigarray(p.wsGenome.chr[id].getnbin(), 0);

  for (auto strand: {Strand::FWD, Strand::REV}) {
    for (auto &x: p.genome.chr[id].getvReadref(strand)) {
      if(x.duplicate) continue;
      addReadToWigArray(p.wsGenome, wigarray, x, p.genome.chr[id].getlen());
    }
  }

  if(p.getMpDir() != "") {
    int binsize(p.wsGenome.getbinsize());
    int mpthre = p.getmpthre() * binsize;
    auto mparray = readMpbl(p.getMpDir(), ("chr" + p.genome.chr[id].getname()), binsize, p.wsGenome.chr[id].getnbin());
    for(int i=0; i<p.wsGenome.chr[id].getnbin(); ++i) {
      if(mparray[i] > mpthre) wigarray.multipleval(i, getratio(binsize, mparray[i]));
    }
  }

  /* Total read normalization */
  if(p.rpm.getType() != "NONE") norm2rpm(p, p.genome.chr[id], wigarray);

  p.wsGenome.setWigStats(id, wigarray);

  /*  t1 = clock();
  clock_t t1,t2;
  p.wsGenome.chr[id].peakcall(wigarray, p.genome.chr[id].getname());
  t2 = clock();
  PrintTime(t1, t2, "peakcall");*/

  return wigarray;
}

void norm2rpm(Mapfile &p, SeqStats &chr, WigArray &wigarray)
{
  static int on(0);
  double w(0);
  std::string ntype(p.rpm.getType());

  if(ntype == "GR") {
    double dn(p.genome.getnread_nonred(Strand::BOTH));
    w = getratio(p.rpm.getnrpm(), dn);
    if(!on) {
      std::cout << boost::format("\ngenomic read number = %1%, after=%2%, w=%3$.3f\n") % (int64_t)dn % p.rpm.getnrpm() % w;
      if(w>2) printwarning(w);
      on=1;
    }
  } else if(ntype == "GD") {
    w = getratio(p.rpm.getndepth(), p.genome.getdepth());
    if(!on) {
      std::cout << boost::format("\ngenomic depth = %1$.2f, after=%2$.2f, w=%3$.3f\n") % p.genome.getdepth() % p.rpm.getndepth() % w;
      if(w>2) printwarning(w);
      on=1;
    }
  } else if(ntype == "CR") {
    double nm = p.rpm.getnrpm() * getratio(chr.getlenmpbl(), p.genome.getlenmpbl());
    double dn = chr.getnread_nonred(Strand::BOTH);
    w = getratio(nm, dn);
    std::cout << boost::format("read number = %1%, after=%2$.1f, w=%3$.3f\n") % static_cast<int64_t>(dn) % nm % w;
    if(w>2) printwarning(w);
  } else if(ntype == "CD") {
    w = getratio(p.rpm.getndepth(), chr.getdepth());
    std::cout << boost::format("depth = %1$.2f, after=%2$.2f, w=%3$.3f\n") % chr.getdepth() % p.rpm.getndepth() % w;
    if(w>2) printwarning(w);
  }

  int32_t nbin(chr.getlen()/p.wsGenome.getbinsize() +1);
  for(int i=0; i<nbin; ++i) { wigarray.multipleval(i, w); }

  chr.setsizefactor(w);
  if(ntype == "GR" || ntype == "GD") p.genome.setsizefactor(w);
  return;
}

void outputWig(Mapfile &p, const std::string &filename)
{
  int32_t binsize(p.wsGenome.getbinsize());

  FILE* File = fopen(filename.c_str(), "w");

  fprintf(File, "track type=wiggle_0\tname=\"%s\"\tdescription=\"Merged tag counts for every %d bp\"\n", p.getSampleName().c_str(), binsize);

  for(size_t i=0; i<p.genome.chr.size(); ++i) {
    WigArray array = makeWigarray(p, i);

    fprintf(File, "variableStep\tchrom=%s\tspan=%d\n", p.genome.chr[i].getrefname().c_str(), binsize);
    array.outputAsWig(File, binsize, p.wsGenome.isoutputzero());
  }
  fclose(File);

  return;
}

void outputBedGraph(Mapfile &p, const std::string &filename)
{
  int32_t binsize(p.wsGenome.getbinsize());

  std::ofstream out(filename);
  out << boost::format("browser position %1%:%2%-%3%\n") % p.genome.chr[1].getrefname() % 0 % (p.genome.chr[1].getlen()/100);
  out << "browser hide all" << std::endl;
  out << "browser pack refGene encodeRegions" << std::endl;
  out << "browser full altGraph" << std::endl;
  out << boost::format("track type=bedGraph name=\"%1%\" description=\"Merged tag counts for every %2% bp\" visibility=full\n")
    % p.getSampleName() % binsize;
  out.close();

  std::string tempfile = filename + ".temp";

  FILE* File = fopen(tempfile.c_str(), "w");

  clock_t t1,t2;
  for(size_t i=0; i<p.genome.chr.size(); ++i) {
    t1 = clock();
    WigArray array = makeWigarray(p, i);
    t2 = clock();
    PrintTime(t1, t2, "makeWigarray");
    t1 = clock();
    array.outputAsBedGraph(File, binsize, p.genome.chr[i].getrefname(), p.genome.chr[i].getlen() -1, p.wsGenome.isoutputzero());
    t2 = clock();
    PrintTime(t1, t2, "outputAsBedGraph");
  }
  fclose (File);

  printf("sort bedGraph...\n");
  std::string command = "sort -k1,1 -k2,2n "+ tempfile +" >> " + filename;
  if(system(command.c_str())) PRINTERR("sorting bedGraph failed.");

  remove(tempfile.c_str());

  return;
}

/*void outputBinary(Mapfile &p, const std::string &filename)
{
  std::ofstream out(filename, std::ios::binary);

  for(size_t i=0; i<p.genome.chr.size(); ++i) {
    WigArray array = makeWigarray(p, i);
    array.outputAsBinary(out);
  }
  return;
}
*/
