/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "pw_gv.h"
#include "pw_readmapfile.h"
#include "util.h"
#include "common.h"
#include "warn.h"
#include "macro.h"
#include "readdata.h"

using namespace std;
using namespace boost::program_options;
namespace fs = boost::filesystem;

variables_map getOpts(int argc, char* argv[]);
void setOpts(options_description &);
void init_dump(const variables_map &);

void printVersion()
{
  cerr << "parse2wig version " << VERSION << endl;
  exit(0);
}  

void help_global() {
    auto helpmsg = R"(
===============

Usage: parse2wig+ [option] -i <inputfile> -o <output> -gt <genome_table>)";

    cerr << "\nparse2wig v" << VERSION << helpmsg << endl;
    return;
}

void checkParam(const variables_map &values){
  vector<string> intopts = {"binsize", "flen", "maxins" , "nrpm", "flen4gc"};
  for (auto x: intopts) chkminus<int>(values, x, 0);
  vector<string> intopts2 = {"rcenter", "thre_pb"};
  for (auto x: intopts2) chkminus<int>(values, x, -1);
  vector<string> dbopts = {"ndepth", "mpthre"};
  for (auto x: dbopts) chkminus<double>(values, x, 0);
  
  if(!RANGE(values["of"].as<int>(), 0, PWFILETYPENUM-1)) printerr("invalid wigfile type.\n");

  string ftype = values["ftype"].as<string>();
  if(ftype != "SAM" && ftype != "BAM" && ftype != "BOWTIE" && ftype != "TAGALIGN") printerr("invalid --ftype.\n");
  string ntype = values["ntype"].as<string>();
  if(ntype != "NONE" && ntype != "GR" && ntype != "GD" && ntype != "CR" && ntype != "CD") printerr("invalid --ntype.\n");

  if (values.count("genome")) {
    if(!values.count("mpbin")) printerr("-GC option requires -mpbin option.\n");
    isFile(values["genome"].as<string>());
  }

  return;
}

variables_map getOpts(int argc, char* argv[])
{
  options_description allopts("Options");
  setOpts(allopts);
  
  variables_map values;
  
  try {
    parsed_options parsed = parse_command_line(argc, argv, allopts);
    store(parsed, values);
    
    if (values.count("version")) printVersion();
    
    if (argc ==1) {
      help_global();
      cerr << "Use --help option for more information on the other options\n\n";
      exit(0);
    }
    if (values.count("help")) {
      help_global();
      cout << "\n" << allopts << endl;
      exit(0);
    }
    vector<string> opts = {"input", "output", "gt"};
    for (auto x: opts) {
      if (!values.count(x)) {
	cerr << "Error: specify --" << x << " option." << endl;
	exit(1);
      }
    }

    notify(values);

    checkParam(values);
    init_dump(values);
  } catch (exception &e) {
    cout << e.what() << endl;
  }
  return values;
}

RefGenome makeRefGenome(variables_map values)
{
  auto gt = read_genometable(values["gt"].as<string>());
  RefGenome g(gt, values["binsize"].as<int>());

  return g;
}

int main(int argc, char* argv[])
{
  variables_map values = getOpts(argc, argv);
  
  fs::path dir(values["odir"].as<string>());
  fs::create_directory(dir);

  RefGenome g = makeRefGenome(values);

#ifdef DEBUG
  g.print();
#endif

    /* 
  if(p->bedfilename){
    isfile(p->bedfilename);
    p->enrichfile = read_bedfile(p->bedfilename, g);
    //    show_bedfile(p->enrichfile, g->chrnum);
  }*/
    
  /* read mapfile */
  read_mapfile(values, g);
  //  Mapfile *mapfile = read_mapfile(p, g);
  
  /* PCR bias filtering and ignore enrichregions */
  /*check_redundant_reads(p, mapfile, g);
  add_read_red_to_genome(mapfile, g);

  if(p->ccp) pw_ccp(p, mapfile, g->chrmax, (int)g->chr[g->chrmax].len);

  if(p->bedfilename) calc_FRiP(p, mapfile, g);*/

  /* calculate depth */
  //  calc_depth(p, mapfile, g);

  /* GC normalization */
  //if(p->genomefile) GCnorm(p, mapfile, g);

  /* make and output wigdata */
  //makewig(p, mapfile, g);

  /* output wigarray_stats, 
     calculate genome coverage */
  // output_wigstats(p, mapfile, g);

  /* output stats */
  // output_stats(p, mapfile, g);

  return 0;
}

void setOpts(options_description &allopts){
  options_description optreq("Required",100);
  optreq.add_options()
    ("input,i",   value<string>(), "Mapping file. Multiple files are allowed (separated by ',')")
    ("output,o",  value<string>(), "Prefix of output files")
    ("gt",        value<string>(), "Genome table (tab-delimited file describing the name and length of each chromosome)")
    ;
  options_description optIO("Input/Output",100);
  optIO.add_options()
    ("binsize,b",   value<int>()->default_value(100),	  "bin size")
    ("ftype,f",     value<string>()->default_value("SAM"), "{SAM|BAM|BOWTIE|TAGALIGN}: format of input file (default:SAM)\nTAGALIGN could be gzip'ed (extension: tagAlign.gz)")
    ("of",        value<int>()->default_value(0),	  "output format\n   0: binary (.bin)\n   1: compressed wig (.wig.gz)\n   2: uncompressed wig (.wig)\n   3: bedGraph (.bedGraph)\n   4: bigWig (.bw)")
    ("odir",        value<string>()->default_value("parse2wigdir"),	  "output directory name")
    ("rcenter", value<int>()->default_value(0), "consider length around the center of fragment ")
    ;
  options_description optsingle("For single-end read",100);
  optsingle.add_options()
    ("flen",        value<int>()->default_value(150), "expected fragment length\n(Automatically calculated in paired-end mode)")
    ;
  options_description optpair("For paired-end read",100);
  optpair.add_options()
    ("pair", 	  "add when the input file is paired-end")
    ("maxins",        value<int>()->default_value(500), "maximum fragment length")
    ;
  options_description optpcr("PCR bias filtering",100);
  optpcr.add_options()
    ("nofilter", 	  "do not filter PCR bias")
    ("thre_pb",        value<int>()->default_value(0),	  "PCRbias threshold (default: more than max(1 read, 10 times greater than genome average)) ")
    ("ncmp",        value<int>()->default_value(10000000),	  "read number for calculating library complexity")
    ;
  options_description optnorm("Total read normalization",100);
  optnorm.add_options()
    ("ntype,n",        value<string>()->default_value("NONE"),  "Total read normalization\n{NONE|GR|GD|CR|CD}\n   NONE: not normalize\n   GR: for whole genome, read number\n   GD: for whole genome, read depth\n   CR: for each chromosome, read number\n   CD: for each chromosome, read depth")
    ("nrpm",        value<int>()->default_value(20000000),	  "Total read number after normalization")
    ("ndepth",      value<double>()->default_value(1.0),	  "Averaged read depth after normalization")
    ("bed",        value<string>(),	  "specify the BED file of enriched regions (e.g., peak regions)")
    ;  
  options_description optmp("Mappability normalization",100);
  optmp.add_options()
    ("mp",        value<string>(),	  "Mappability file")  // mpbin ni touitu?
    ("mpthre",    value<double>()->default_value(0.3),	  "Threshold of low mappability regions")
    ("mpbin",        value<string>(),	  "mpbinaryfile")
    ;
  options_description optgc("GC bias normalization\n   (require large time and memory)",100);
  optgc.add_options()
    ("genome",     value<string>(),	  "reference genome sequence for GC content estimation")
    ("flen4gc",    value<int>()->default_value(120),  "fragment length for calculation of GC distribution")
    ("gcdepthoff", "do not consider depth of GC contents")
    ;
  options_description optother("Others",100);
  optother.add_options()
    ("ccp", 	  "make cross-correlation profile")
    ("version,v", "print version")
    ("help,h", "show help message")
    ;
  allopts.add(optreq).add(optIO).add(optsingle).add(optpair).add(optpcr).add(optnorm).add(optmp).add(optgc).add(optother);
  return;
}

void init_dump(const variables_map &values){
  vector<string> str_wigfiletype = {"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};

  BPRINT("\n======================================\n");
  BPRINT("parse2wig version %1%\n\n") % VERSION;
  BPRINT("Input file %1%\n")         % values["input"].as<string>();
  BPRINT("\tFormat: %1%\n")          % values["ftype"].as<string>();
  BPRINT("Output file: %1%/%2%\n")   % values["odir"].as<string>() % values["output"].as<string>();
  BPRINT("\tFormat: %1%\n")          % str_wigfiletype[values["of"].as<int>()];
  BPRINT("Genome-table file: %1%\n") % values["gt"].as<string>();
  BPRINT("Binsize: %1% bp\n")        % values["binsize"].as<int>();
  if (!values.count("pair")) BPRINT("Predefined fragment length: %1%\n") % values["flen"].as<int>();
  else BPRINT("Maximum fragment length: %1%\n") % values["maxins"].as<int>();
  if (!values.count("nofilter")) {
    BPRINT("PCR bias filtering: ON\n");
    if (values["thre_pb"].as<int>()) BPRINT("PCR bias threshold: > %1%\n") % values["thre_pb"].as<int>();
  }
  BPRINT("\t%1% reads used for library complexity\n") % values["ncmp"].as<int>();
  if (values.count("bed")) BPRINT("Bed file: %1%\n") % values["bed"].as<string>();
  if (values.count("ccp")) BPRINT("Make cross-correlation profile.\n");

  string ntype = values["ntype"].as<string>();
  BPRINT("\nTotal read normalization: %1%\n") % ntype;
  if(ntype == "GR" || ntype == "CR"){
    BPRINT("\tnormed read: %1% M for genome\n") % (values["nrpm"].as<int>() /(double)NUM_1M);
  }
  else if(ntype == "GD" || ntype == "CD"){
    BPRINT("\tnormed depth: %1%\n") % values["ndepth"].as<double>();
  }
  printf("\n");
  if (values.count("mp")) {
    printf("Mappability normalization:\n");
    BPRINT("\tfile prefix: %1%\n") % values["mp"].as<string>();
    BPRINT("\tLow mappablitiy threshold: %1%\n") % values["mpthre"].as<double>();
  }
  if (values.count("genome")) {
    printf("Correcting GC bias:\n");
    BPRINT("\tChromosome directory: %1%\n") % values["genome"].as<string>();
    BPRINT("\tmappability binary prefix: %1%\n") % values["mpbin"].as<string>();
    BPRINT("\tLength for GC distribution: %1%\n") % values["flen4gc"].as<int>();
  }
  printf("======================================\n");
  return;
}
