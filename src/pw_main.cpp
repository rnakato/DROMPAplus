/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <boost/filesystem.hpp>
#include "SSP/src/ParseMapfile.hpp"
#include "pw_makefile.h"
#include "pw_gc.h"
#include "version.h"
#include "pw_gv.h"
#include "SSP/src/ssp_shiftprofile.h"
#include "readbpstatus.h"
#include "mytype.h"

namespace {
  const int numGcov(5000000);
}

MyOpt::Variables getOpts(Mapfile &p, int argc, char* argv[]);
void setOpts(MyOpt::Opts &);
void init_dump(const MyOpt::Variables &);
void output_stats(const MyOpt::Variables &values, const Mapfile &p);
void calcGenomeCoverage(const MyOpt::Variables &values, Mapfile &p);
void output_wigstats(Mapfile &p);

void printVersion()
{
  std::cerr << "parse2wig version " << VERSION << std::endl;
  exit(0);
}

void help_global()
{
  auto helpmsg = R"(
===============

Usage: parse2wig+ [option] -i <inputfile> -o <output> -gt <genome_table>)";
  
  std::cerr << "\nparse2wig v" << VERSION << helpmsg << std::endl;
  return;
}

void calcFRiP(SeqStats &chr, const std::vector<bed> &vbed)
{
  uint64_t nread_inbed(0);
  std::vector<BpStatus> array(chr.getlen(), BpStatus::MAPPABLE);
  arraySetBed(array, chr.getname(), vbed);
  for (auto strand: {Strand::FWD, Strand::REV}) {
    for (auto &x: chr.getvReadref_notconst(strand)) {
      if(x.duplicate) continue;
      int s(std::min(x.F3, x.F5));
      int e(std::max(x.F3, x.F5));
      for(int i=s; i<=e; ++i) {
	if(array[i] == BpStatus::INBED) {
	  x.inpeak = 1;
	  ++nread_inbed;
	  break;
	}
      }
    }
  }
  chr.setFRiP(nread_inbed);
  return;
}

void setFRiP(SeqStatsGenome &genome)
{
  std::cout << "calculate FRiP score.." << std::flush;
  for(auto &x: genome.chr) calcFRiP(x, genome.getvbedref());
  std::cout << "done." << std::endl;
  return;
}

template <class T>
void calcdepth(T &obj, const int32_t flen)
{
  uint64_t lenmpbl = obj.getlenmpbl();
  double d = lenmpbl ? getratio(obj.getnread_nonred(Strand::BOTH) * flen, lenmpbl): 0;
  obj.setdepth(d);
}

int main(int argc, char* argv[])
{
  Mapfile p;
  MyOpt::Variables values = getOpts(p, argc, argv);
  p.setValues(values);

  read_mapfile(p.genome);

  p.genome.dflen.outputDistFile(p.getprefix(), p.genome.getnread(Strand::BOTH));

  p.complexity.checkRedundantReads(p.genome);
 
  strShiftProfile(p.sspst, p.genome, p.getprefix(), "jaccard");
  for (auto &x: p.genome.chr) calcdepth(x, p.genome.dflen.getflen());
  calcdepth(p.genome, p.genome.dflen.getflen());
  

  // BED file
  if (values.count("bed")) {
    p.genome.setbed(values["bed"].as<std::string>());
    setFRiP(p.genome);
  }

#ifdef DEBUG
  p.genome.printReadstats();
#endif

  // Genome coverage
  calcGenomeCoverage(values, p);
  
  // GC contents
  if (values.count("genome")) normalizeByGCcontents(values, p);
  
  // make and output wigdata
  makewig(values, p);
  
  p.estimateZINB();  // for genome
  
  p.printPeak(values);
  
  // output stats
  output_wigstats(p);
  output_stats(values, p);

  return 0;
}

void checkParam(const MyOpt::Variables &values)
{
  std::vector<std::string> intopts = {"binsize", "flen", "maxins" , "nrpm", "flen4gc", "threads"};
  for (auto x: intopts) chkminus<int>(values, x, 0);
  std::vector<std::string> intopts2 = {"rcenter", "thre_pb"};
  for (auto x: intopts2) chkminus<int>(values, x, -1);
  std::vector<std::string> dbopts = {"ndepth", "mpthre"};
  for (auto x: dbopts) chkminus<double>(values, x, 0);
  
  if(!my_range(values["of"].as<int>(), 0, static_cast<int>(WigType::WIGTYPENUM) -1)) PRINTERR("invalid wigfile type.\n");

  if(values.count("ftype")) {
    std::string ftype = values["ftype"].as<std::string>();
    if(ftype != "SAM" && ftype != "BAM" && ftype != "BOWTIE" && ftype != "TAGALIGN") PRINTERR("invalid --ftype.\n");
  }
  std::string ntype = values["ntype"].as<std::string>();
  if(ntype != "NONE" && ntype != "GR" && ntype != "GD" && ntype != "CR" && ntype != "CD") PRINTERR("invalid --ntype.\n");

  if (values.count("genome")) {
    //    if(!values.count("mp")) PRINTERR("--genome option requires --mp option.\n");
    isFile(values["genome"].as<std::string>());
  }
  
  return;
}

MyOpt::Variables getOpts(Mapfile &p, int argc, char* argv[])
{
  DEBUGprint("setOpts...");

  MyOpt::Opts allopts("Options");
  p.setOpts(allopts);
  setOpts(allopts);
  
  MyOpt::Variables values;
  
  DEBUGprint("getOpts...");

  try {
    boost::program_options::parsed_options parsed = parse_command_line(argc, argv, allopts);
    store(parsed, values);
    
    if (values.count("version")) printVersion();

    if (argc ==1) {
      help_global();
      std::cerr << "Use --help option for more information on the other options\n\n";
      exit(0);
    }
    if (values.count("help")) {
      help_global();
      std::cout << "\n" << allopts << std::endl;
      exit(0);
    }
    std::vector<std::string> opts = {"input", "output", "gt"};
    for (auto x: opts) {
      if (!values.count(x)) PRINTERR("specify --" << x << " option.");
    }

    notify(values);
    checkParam(values);

    boost::filesystem::path dir(values["odir"].as<std::string>());
    boost::filesystem::create_directory(dir);

    init_dump(values);
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    exit(0);
  }
  
  DEBUGprint("getOpts done.");
  return values;
}

void setOpts(MyOpt::Opts &allopts)
{
  using namespace boost::program_options;

  MyOpt::setOptIO(allopts, "parse2wigdir+");
  MyOpt::setOptPair(allopts);
  MyOpt::Opts optOut("Output",100);
  optOut.add_options()
    ("of",        value<int>()->default_value(0), "output format\n   0: binary (.bin)\n   1: compressed wig (.wig.gz)\n   2: uncompressed wig (.wig)\n   3: bedGraph (.bedGraph)\n   4: bigWig (.bw)")
    ("rcenter", value<int>()->default_value(0), "consider length around the center of fragment ")
    ;
  MyOpt::Opts optnorm("Total read normalization",100);
  optnorm.add_options()
    ("ntype,n",        value<std::string>()->default_value("NONE"),  "Total read normalization\n{NONE|GR|GD|CR|CD}\n   NONE: not normalize\n   GR: for whole genome, read number\n   GD: for whole genome, read depth\n   CR: for each chromosome, read number\n   CD: for each chromosome, read depth")
    ("nrpm",        value<int>()->default_value(2*NUM_10M),	  "Total read number after normalization")
    ("ndepth",      value<double>()->default_value(1.0),	  "Averaged read depth after normalization")
    ("bed",        value<std::string>(),	  "specify the BED file of enriched regions (e.g., peak regions)")
    ;  
  MyOpt::Opts optmp("Mappability normalization",100);
  optmp.add_options()
    ("mp",        value<std::string>(),	  "Mappability file")
    ("mpthre",    value<double>()->default_value(0.3),	  "Threshold of low mappability regions")
    ;
  MyOpt::Opts optgc("GC bias normalization\n   (require large time and memory)",100);
  optgc.add_options()
    ("genome",     value<std::string>(),	  "reference genome sequence for GC content estimation")
    ("flen4gc",    value<int>()->default_value(120),  "fragment length for calculation of GC distribution")
    ("gcdepthoff", "do not consider depth of GC contents")
    ;
  allopts.add(optOut).add(optnorm).add(optmp).add(optgc);
  MyOpt::setOptOther(allopts);
  return;
}

void init_dump(const MyOpt::Variables &values){
  std::vector<std::string> str_wigfiletype = {"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};
 
  std::cout << boost::format("\n======================================\n");
  std::cout << boost::format("parse2wig version %1%\n\n") % VERSION;

  MyOpt::dumpIO(values);
  MyOpt::dumpGenomeTable(values);

  std::cout << boost::format("\tOutput format: %1%\n")          % str_wigfiletype[values["of"].as<int>()];
  std::cout << boost::format("Binsize: %1% bp\n")        % values["binsize"].as<int>();

  MyOpt::dumpFragmentLengthDist(values);
  MyOpt::dumpPair(values);
  MyOpt::dumpLibComp(values);

  if (values.count("bed")) std::cout << boost::format("Bed file: %1%\n")  % values["bed"].as<std::string>();

  std::string ntype = values["ntype"].as<std::string>();
  std::cout << boost::format("\nTotal read normalization: %1%\n") % ntype;
  if(ntype == "GR" || ntype == "CR"){
    std::cout << boost::format("\tnormed read: %1% M for genome\n") % getratio(values["nrpm"].as<int>(), NUM_1M);
  }
  else if(ntype == "GD" || ntype == "CD"){
    std::cout << boost::format("\tnormed depth: %1%\n") % values["ndepth"].as<double>();
  }
  printf("\n");
  if (values.count("mp")) {
    printf("Mappability normalization:\n");
    std::cout << boost::format("\tFile directory: %1%\n") % values["mp"].as<std::string>();
    std::cout << boost::format("\tLow mappablitiy threshold: %1%\n") % values["mpthre"].as<double>();
  }
  if (values.count("genome")) {
    printf("Correcting GC bias:\n");
    std::cout << boost::format("\tChromosome directory: %1%\n") % values["genome"].as<std::string>();
  }
  MyOpt::dumpOther(values);
  printf("======================================\n");
  return;
}

template <class T>
void print_SeqStats(const MyOpt::Variables &values, std::ofstream &out, const T &p, const Mapfile &mapfile)
{
  /* genome data */
  out << p.getname() << "\t" << p.getlen()  << "\t" << p.getlenmpbl() << "\t" << p.getpmpbl() << "\t";
  /* total reads*/
  out << boost::format("%1%\t%2%\t%3%\t%4$.1f%%\t")
    % p.getnread(Strand::BOTH) % p.getnread(Strand::FWD) % p.getnread(Strand::REV)
    % getpercent(p.getnread(Strand::BOTH), mapfile.genome.getnread(Strand::BOTH));

  std::vector<Strand::Strand> vstr = {Strand::BOTH, Strand::FWD, Strand::REV};
  for (auto strand: vstr) printNumandPer(out, p.getnread_nonred(strand), p.getnread(strand));
  for (auto strand: vstr) printNumandPer(out, p.getnread_red(strand),    p.getnread(strand));

  /* reads after GCnorm */
  if(values.count("genome")) {
    for (auto strand: vstr) printNumandPer(out, p.getnread_afterGC(strand), p.getnread(strand));
  }
  out << boost::format("%1$.3f\t") % p.getdepth();
  if(p.getsizefactor()) out << boost::format("%1$.3f\t") % p.getsizefactor();
  else                  out << " - \t";
  if(values["ntype"].as<std::string>() == "NONE") out << p.getnread_nonred(Strand::BOTH) << "\t";
  else out << p.getnread_rpm(Strand::BOTH) << "\t";
  
  if(mapfile.islackOfRead4GenomeCov()) {
    out << boost::format("%1$.3f\t(%2$.3f)\t")
      % getratio(p.getncov(),     p.getnbp())
      % getratio(p.getncovnorm(), p.getnbp());
  } else {
    out << boost::format("%1$.3f\t%2$.3f\t")
      % getratio(p.getncov(),     p.getnbp())
      % getratio(p.getncovnorm(), p.getnbp());
  }

  if(values.count("bed")) out << boost::format("%1$.3f\t") % p.getFRiP();

  return;
}

void output_stats(const MyOpt::Variables &values, const Mapfile &p)
{
  std::string filename = p.getbinprefix() + ".csv";
  std::ofstream out(filename);

  out << "parse2wig version " << VERSION << std::endl;
  out << "Input file: \"" << values["input"].as<std::string>() << "\"" << std::endl;
  out << "Redundancy threshold: >" << p.complexity.getThreshold() << std::endl;

  p.complexity.print(out);
  p.genome.dflen.printFlen(out);
  if(values.count("genome")) out << "GC summit: " << p.getmaxGC() << std::endl;

  // Global stats
  out << "\n\tlength\tmappable base\tmappability\t";
  out << "total reads\t\t\t\t";
  out << "nonredundant reads\t\t\t";
  out << "redundant reads\t\t\t";
  if(values.count("genome")) out << "reads (GCnormed)\t\t\t";
  out << "read depth\t";
  out << "scaling weight\t";
  out << "normalized read number\t";
  out << "gcov (Raw)\tgcov (Normed)\t";
  out << "bin mean\tbin variance\t";
  if(values.count("bed")) out << "FRiP\t";
  out << "nb_p\tnb_n\tnb_p0\t";
  out << std::endl;
  out << "\t\t\t\t";
  out << "both\tforward\treverse\t% genome\t";
  out << "both\tforward\treverse\t";
  out << "both\tforward\treverse\t";
  if(values.count("genome")) out << "both\tforward\treverse\t";
  out << std::endl;

  // SeqStats
  print_SeqStats(values, out, p.genome, p);
  p.wsGenome.printPoispar(out);
  p.wsGenome.printZINBpar(out);
  out << std::endl;
  
  for(auto &chr: p.genome.chr) {
    print_SeqStats(values, out, chr, p);
    chr.ws.printPoispar(out);
    chr.ws.printZINBpar(out);
    out << std::endl;
  }
  
  std::cout << "stats is output in " << filename << "." << std::endl;

  return;
}

std::vector<BpStatus> makeGcovArray(const MyOpt::Variables &values, SeqStats &chr, Mapfile &p, double r4cmp)
{
  std::vector<BpStatus> array;
  if(values.count("mp")) array = readMpbl_binary(values["mp"].as<std::string>(), ("chr" + p.lchr->getname()), chr.getlen());
  else array = readMpbl_binary(chr.getlen());
  if(values.count("bed")) arraySetBed(array, chr.getname(), p.genome.getvbedref());

  int32_t size = array.size();
  for (auto strand: {Strand::FWD, Strand::REV}) {
    const std::vector<Read> &vReadref = chr.getvReadref(strand);
    for (auto &x: vReadref) {
      if(x.duplicate) continue;
      
      BpStatus val(BpStatus::UNMAPPABLE);
      if(rand() >= r4cmp) val = BpStatus::COVREAD_ALL; else val = BpStatus::COVREAD_NORM;
      
      int32_t s(std::max(0, std::min(x.F3, x.F5)));
      int32_t e(std::min(std::max(x.F3, x.F5), size-1));
      if(s >= size || e < 0) {
	std::cerr << "Warning: " << chr.getname() << " read " << s <<"-"<< e << " > array size " << array.size() << std::endl;
      }
      for(int32_t i=s; i<=e; ++i) if(array[i]==BpStatus::MAPPABLE) array[i]=val;
    }
  }
  return array;
}

void calcGcovchr(const MyOpt::Variables &values, Mapfile &p, int32_t s, int32_t e, double r4cmp)
{
  for(int32_t i=s; i<=e; ++i) {
    std::cout << p.genome.chr[i].getname() << ".." << std::flush;
    auto array = makeGcovArray(values, p.genome.chr[i], p, r4cmp);
    p.genome.chr[i].calcGcov(array);
  }
}

void calcGenomeCoverage(const MyOpt::Variables &values, Mapfile &p)
{
  std::cout << "calculate genome coverage.." << std::flush;

  // ignore peak region
  double r = getratio(numGcov, p.genome.getnread_nonred(Strand::BOTH) - p.genome.getnread_inbed());
  if(r>1){
    std::cerr << "Warning: number of reads is < "<< static_cast<int>(numGcov/NUM_1M) << " million.\n";
    p.lackOfRead4GenomeCov_on();
  }
  double r4cmp = r*RAND_MAX;

  boost::thread_group agroup;
  for(uint i=0; i<p.genome.vsepchr.size(); i++) {
    agroup.create_thread(bind(calcGcovchr, boost::cref(values), boost::ref(p), p.genome.vsepchr[i].s, p.genome.vsepchr[i].e, r4cmp));
  }
  agroup.join_all();
  
  std::cout << "done." << std::endl;
  return;
}

void output_wigstats(Mapfile &p)
{
  std::string filename = p.getbinprefix() + ".binarray_dist.csv";
  std::ofstream out(filename);

  std::cout << "generate " << filename << ".." << std::flush;

  out << "\tGenome\t\t\t";
  for (auto &x: p.genome.chr) out << x.getname() << "\t\t\t\t";
  out << std::endl;
  out << "read number\tnum of bins genome\tprop\tZINB estimated\t";
  for (size_t i=0; i<p.genome.chr.size(); ++i) out << "num of bins\tprop\tPoisson estimated\tZINB estimated\t";
  out << std::endl;

  for(size_t i=0; i<p.wsGenome.wigDist.size(); ++i) {
    out << i << "\t";
    p.wsGenome.printwigDist(out, i);
    out << p.wsGenome.getZINB(i) << "\t";
    //    out << p.genome.getZIP(i) << "\t";
    for (auto &x: p.genome.chr) {
      x.ws.printwigDist(out, i);
      out << x.ws.getPoisson(i) << "\t";
      out << x.ws.getZINB(i) << "\t";
    }
    out << std::endl;
  }

  std::cout << "done." << std::endl;
  return;
}
