/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "pw_makefile.hpp"
#include "pw_gv.hpp"
#include "WigStats.hpp"
#include "ReadMpbldata.hpp"
#include "../submodules/SSP/src/SeqStats.hpp"

namespace {
  void printwarning(double w)
  {
    std::cerr << "Warning: Read scaling weight = "
              << w
              << ". Too much scaling up will bring noisy results."
              << std::endl;
  }

  void addReadToWigArray(const WigStatsGenome &p, WigArray &wigarray, const Read x, const int64_t chrlen, const int32_t readlenF3, const int32_t readlenF5)
  {
    int32_t s, e;
    s = std::min(x.F3, x.F5);
    e = std::max(x.F3, x.F5);

    int32_t rcenter(p.getrcenter());
    if (rcenter) {  // consider only center region of fragments
      s = (s + e - rcenter)/2;
      e = s + rcenter;
    }
    s = std::max(0, s);
    e = std::min(e, (int32_t)(chrlen -1));

    if (p.isonlyreadregion() && (e-s) > 300) { // for paired-end: consider only read region
      int32_t sbin(s/p.getbinsize());
      int32_t ebin((e+readlenF3)/p.getbinsize());
      for (int32_t j=sbin; j<=ebin; ++j) wigarray.addval(j, x.getWeight());
      sbin = (e-readlenF5)/p.getbinsize();
      ebin = e/p.getbinsize();
      for (int32_t j=sbin; j<=ebin; ++j) wigarray.addval(j, x.getWeight());
    } else {
      int32_t sbin(s/p.getbinsize());
      int32_t ebin(e/p.getbinsize());
      for (int32_t j=sbin; j<=ebin; ++j) wigarray.addval(j, x.getWeight());
    }
    return;
  }

  double getScaleWeight_for_totalreads(Mapfile &p, const SeqStats &chr)
  {
    static int32_t on(0);
    double w(0);
    std::string ntype(p.rpm.getType());

    if (ntype == "GR") {
      double dn(p.genome.getnread_nonred(Strand::BOTH));
      w = getratio(p.rpm.getnrpm(), dn);
      if (!on) {
        std::cout << boost::format("\ngenomic read number = %1%, after=%2%, w=%3$.3f\n") % (int64_t)dn % p.rpm.getnrpm() % w;
        if (w>2) printwarning(w);
        on=1;
      }
    } else if (ntype == "GD") {
      w = getratio(p.rpm.getndepth(), p.genome.getdepth());
      if (!on) {
        std::cout << boost::format("\ngenomic depth = %1$.2f, after=%2$.2f, w=%3$.3f\n") % p.genome.getdepth() % p.rpm.getndepth() % w;
        if (w>2) printwarning(w);
        on=1;
      }
    } else if (ntype == "CR") {
      double nm = p.rpm.getnrpm() * getratio(chr.getlenmpbl(), p.genome.getlenmpbl());
      double dn = chr.getnread_nonred(Strand::BOTH);
      w = getratio(nm, dn);
      std::cout << boost::format("read number = %1%, after=%2$.1f, w=%3$.3f\n") % static_cast<int64_t>(dn) % nm % w;
      if (w>2) printwarning(w);
    } else if (ntype == "CD") {
      w = getratio(p.rpm.getndepth(), chr.getdepth());
      std::cout << boost::format("depth = %1$.2f, after=%2$.2f, w=%3$.3f\n") % chr.getdepth() % p.rpm.getndepth() % w;
      if (w>2) printwarning(w);
    }

    return w;
  }

  WigArray count_and_normalize_Wigarray(Mapfile &p, const int32_t id)
  {
    std::cout << "chr" << p.genome.chr[id].getname() << ".." << std::flush;
    WigArray wigarray(p.wsGenome.chr[id].getnbin(), 0);

    // Convert readarray to Wig
    for (auto strand: {Strand::FWD, Strand::REV}) {
      for (auto &x: p.genome.chr[id].getvReadref(strand)) {
        if (x.duplicate) continue;
        addReadToWigArray(p.wsGenome, wigarray, x, p.genome.chr[id].getlen(), p.genome.dflen.getlenF3(), p.genome.dflen.getlenF5());
      }
    }

    // Mappability normalization
    if (p.getMpblBinaryDir() != "") {
      int32_t binsize(p.wsGenome.getbinsize());
      int32_t mpthre = p.getmpthre() * binsize;
      auto mparray = readMpblWigArray(p.getMpblBinaryDir(),
                                      ("chr" + p.genome.chr[id].getname()),
                                      binsize,
                                      p.wsGenome.chr[id].getnbin());
      for (int32_t i=0; i<p.wsGenome.chr[id].getnbin(); ++i) {
        //      std::cout << "mparray[i]: " << mparray[i] << std::endl;
        if (mparray[i] > mpthre) wigarray.multipleval(i, getratio(binsize, mparray[i]));
      }
    }

    /* Total read normalization */
    if (p.rpm.getType() != "NONE") {
      double w = getScaleWeight_for_totalreads(p, p.genome.chr[id]);
      p.genome.setsizefactor(w, id);
      if (p.rpm.getType() == "GR" || p.rpm.getType() == "GD") p.genome.setsizefactor(w);

      for (int32_t i=0; i<p.wsGenome.chr[id].getnbin(); ++i) { wigarray.multipleval(i, w); }
    }

    p.wsGenome.setWigStats(id, wigarray);

    // Peak calling
    /*  t1 = clock();
        clock_t t1,t2;
        p.wsGenome.chr[id].peakcall(wigarray, p.genome.chr[id].getname());
        t2 = clock();
        PrintTime(t1, t2, "peakcall");*/

    return wigarray;
  }

  void outputWig(Mapfile &p, const std::string &filename)
  {
    int32_t binsize(p.wsGenome.getbinsize());

    FILE* File = fopen(filename.c_str(), "w");

    fprintf(File, "track type=wiggle_0\tname=\"%s\"\tdescription=\"Merged tag counts for every %d bp\"\n", p.getSampleName().c_str(), binsize);

    for (size_t i=0; i<p.genome.getnchr(); ++i) {
      WigArray array = count_and_normalize_Wigarray(p, i);

      fprintf(File, "variableStep\tchrom=%s\tspan=%d\n", p.genome.chr[i].getrefname().c_str(), binsize);
      bool isfloat(false);
      array.outputAsWig(File, binsize, p.wsGenome.isoutputzero(), isfloat);
    }
    fclose(File);

    return;
  }

  void outputBedGraph(Mapfile &p, const std::string &filename)
  {
    int32_t binsize(p.wsGenome.getbinsize());
  printf("ss array: \n");
    std::ofstream out(filename);
    out << boost::format("browser position %1%:%2%-%3%\n") % p.genome.chr[0].getrefname() % 0 % (p.genome.chr[0].getlen()/100);
    out << "browser hide all" << std::endl;
    out << "browser pack refGene encodeRegions" << std::endl;
    out << "browser full altGraph" << std::endl;
    out << boost::format("track type=bedGraph name=\"%1%\" description=\"Merged tag counts for every %2% bp\" visibility=full\n")
      % p.getSampleName() % binsize;
    out.close();

  printf("ss afffrray: \n");
    std::string tempfile = filename + ".tmpfile";

    FILE* File = fopen(tempfile.c_str(), "w");

    clock_t t1,t2;
    for (size_t i=0; i<p.genome.getnchr(); ++i) {
      t1 = clock();
      WigArray array = count_and_normalize_Wigarray(p, i);
      t2 = clock();
      PrintTime(t1, t2, "count_and_normalize_Wigarray");
      t1 = clock();
      bool isfloat(false);
      array.outputAsBedGraph(File,
                             binsize,
                             p.genome.chr[i].getrefname(),
                             p.genome.chr[i].getlen() -1,
                             p.wsGenome.isoutputzero(),
                             isfloat);
      t2 = clock();
      PrintTime(t1, t2, "outputAsBedGraph");
    }
    fclose (File);

    printf("sort bedGraph...\n");
    std::string command = "sort -k1,1 -k2,2n "+ tempfile +" >> " + filename;
    if (system(command.c_str())) PRINTERR_AND_EXIT("sorting bedGraph failed.");

    remove(tempfile.c_str());

    return;
  }

}
void generate_wigfile(Mapfile &p)
{
  printf("Convert read data to array: \n");
  WigType oftype(p.wsGenome.getWigType());
  std::string filename(p.getbinprefix());

  if (oftype==WigType::COMPRESSWIG || oftype==WigType::UNCOMPRESSWIG) {
    filename += ".wig";
    outputWig(p, filename);
    if (oftype==WigType::COMPRESSWIG) {
      std::string command = "gzip -f " + filename;
      if (system(command.c_str())) PRINTERR_AND_EXIT("gzip .wig failed.");
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
    outputBedGraph(p, std::string(tmpfile));
    printf("Convert to bigWig...\n");
    std::string command = "bedGraphToBigWig " + std::string(tmpfile) + " " + p.genome.getGenomeTable() + " " + p.getbinprefix() + ".bw";
    if (system(command.c_str())) {
      unlink(tmpfile);
      std::cerr << "Error: command " << command << "return nonzero status. "
                << "Add the PATH to 'DROMPAplus/otherbins'." << std::endl;
    }
    unlink(tmpfile);
  }

  printf("done.\n");
  return;
}
