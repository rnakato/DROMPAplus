/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "pw_readmapfile.h"
#include "readdata.h"
#include "pw_makefile.h"
#include "pw_gv.h"
#include "pw_gc.h"

variables_map getOpts(int argc, char* argv[]);
void setOpts(options_description &);
void init_dump(const variables_map &);
void output_stats(const variables_map &values, const Mapfile &p);
void calcGenomeCoverage(const variables_map &values, Mapfile &p);
void output_wigstats(const variables_map &values, Mapfile &p);

Mapfile::Mapfile(const variables_map &values):
  yeast(false), genome("Genome"), flen_def(values["flen"].as<int>()), thre4filtering(0), nt_all(0), nt_nonred(0), nt_red(0), tv(0), gv(0), r4cmp(0), maxGC(0)
{
  oprefix = values["odir"].as<string>() + "/" + values["output"].as<string>();

  readGenomeTable(values);
  if(values.count("mp")) getMpbl(values["mp"].as<string>(), chr);

  for(auto &x:chr) {
    x.nbin = x.getlen()/values["binsize"].as<int>() +1;
    genome.addlen(x);
  }

  // yeast
  for(auto x:chr) if(x.name == "I" || x.name == "II" || x.name == "III") yeast = true;
  for(auto &x:chr) if(yeast) x.yeaston();

  vector<double> gcw(values["flen4gc"].as<int>(),0);
  GCweight = gcw;

  // sepchr
  vsepchr = getVsepchr(values["threads"].as<int>());

#ifdef DEBUG
  cout << "chr\tautosome" << endl;
  for(auto x:chr) {
    cout << x.name << "\t" << x.isautosome() << endl;
  }
  for(uint i=0; i<vsepchr.size(); i++) cout << "thread " << (i+1) << ": "<< vsepchr[i].s << "-" << vsepchr[i].e << endl;

  printstats();

#endif
}

void Mapfile::readGenomeTable(const variables_map &values)
{
  vector<string> v;
  string lineStr;
  ifstream in(values["gt"].as<string>());
  if(!in) PRINTERR("Could nome open " << values["gt"].as<string>() << ".");

  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0] == '#') continue;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    SeqStats s(v[0], stoi(v[1]));    
    chr.push_back(s);
  }
  
  long lenmax(0);
  for(auto itr = chr.begin(); itr != chr.end(); ++itr) {
    if(lenmax < itr->getlen()) {
      lenmax = itr->getlen();
      lchr = itr;
    }
  }
  return;
}

vector<sepchr> Mapfile::getVsepchr(const int numthreads)
{
  vector<sepchr> vsepchr;

  uint sepsize = genome.getlen()/numthreads;
  for(uint i=0; i<chr.size(); ++i) {
    uint s = i;
    long len(0);
    while(len < sepsize && i<chr.size()) {
      len += chr[i].getlen();
      i++;
    }
    i--;
    uint e = i;
    sepchr sep(s,e);
    vsepchr.push_back(sep);
  }
  return vsepchr;
}

void getMpbl(const string mpdir, vector<SeqStats> &chr)
{
  string lineStr;
  vector<string> v;
  string mpfile = mpdir + "/map_fragL150_genome.txt";
  ifstream in(mpfile);
  if(!in) PRINTERR("Could nome open " << mpfile << ".");
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0] == '#') continue;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    for(auto &x: chr) {
      //     setlenmpbl(x, rmchr(v[0]), stoi(v[1]));
      if(x.name == rmchr(v[0])) x.len_mpbl = stoi(v[1]);
    }
  }
  return;
}

void printVersion()
{
  cerr << "parse2wig version " << VERSION << endl;
  exit(0);
}

void help_global()
{
  auto helpmsg = R"(
===============

Usage: parse2wig+ [option] -i <inputfile> -o <output> -gt <genome_table>)";
  
  cerr << "\nparse2wig v" << VERSION << helpmsg << endl;
  return;
}

void outputPeak(const variables_map &values, Mapfile &p)
{
  string filename = p.oprefix + "." + IntToString(values["binsize"].as<int>()) + ".peak.xls";
  int binsize = values["binsize"].as<int>();
  ofstream out(filename);

  p.vPeak[0].printHead(out);
  for(uint i=0; i<p.vPeak.size(); ++i) {
    p.vPeak[i].print(out, i, binsize);
  }
  return;
}

int main(int argc, char* argv[])
{
  variables_map values = getOpts(argc, argv);
  
  boost::filesystem::path dir(values["odir"].as<string>());
  boost::filesystem::create_directory(dir);

  Mapfile p(values);
  read_mapfile(values, p);

  // PCR bias filtering
  check_redundant_reads(values, p);
  p.setnread_red();

  estimateFragLength(values, p);

  // BED file
  if (values.count("bed")) {
    string bedfile = values["bed"].as<string>();
    isFile(bedfile);
    p.vbed = parseBed<bed>(bedfile);
    //    printBed(p.vbed);
    p.setFRiP();
  }

#ifdef DEBUG
  p.printstats();
#endif

  // Genome coverage
  calcGenomeCoverage(values, p);
  
  // GC contents
  if (values.count("genome")) {
    make_GCdist(values, p);
    weightRead(values, p);
  }

  // make and output wigdata
  makewig(values, p);
  p.estimateZINB();  // for genome
  
  outputPeak(values, p);
  
  // output stats
  output_wigstats(values, p);
  output_stats(values, p);

  return 0;
}

void checkParam(const variables_map &values)
{
  vector<string> intopts = {"binsize", "flen", "maxins" , "nrpm", "flen4gc", "threads"};
  for (auto x: intopts) chkminus<int>(values, x, 0);
  vector<string> intopts2 = {"rcenter", "thre_pb"};
  for (auto x: intopts2) chkminus<int>(values, x, -1);
  vector<string> dbopts = {"ndepth", "mpthre"};
  for (auto x: dbopts) chkminus<double>(values, x, 0);
  
  if(!RANGE(values["of"].as<int>(), 0, PWFILETYPENUM-1)) PRINTERR("invalid wigfile type.\n");

  string ftype = values["ftype"].as<string>();
  if(ftype != "SAM" && ftype != "BAM" && ftype != "BOWTIE" && ftype != "TAGALIGN") PRINTERR("invalid --ftype.\n");
  string ntype = values["ntype"].as<string>();
  if(ntype != "NONE" && ntype != "GR" && ntype != "GD" && ntype != "CR" && ntype != "CD") PRINTERR("invalid --ntype.\n");

  if (values.count("genome")) {
    //    if(!values.count("mp")) PRINTERR("--genome option requires --mp option.\n");
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
      if (!values.count(x)) PRINTERR("specify --" << x << " option.");
    }

    notify(values);

    checkParam(values);
    init_dump(values);
  } catch (exception &e) {
    cout << e.what() << endl;
    exit(0);
  }
  return values;
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
    ("binsize,b",   value<int>()->default_value(50),	  "bin size")
    ("ftype,f",     value<string>()->default_value("SAM"), "{SAM|BAM|BOWTIE|TAGALIGN}: format of input file (default:SAM)\nTAGALIGN could be gzip'ed (extension: tagAlign.gz)")
    ("of",        value<int>()->default_value(0),	  "output format\n   0: binary (.bin)\n   1: compressed wig (.wig.gz)\n   2: uncompressed wig (.wig)\n   3: bedGraph (.bedGraph)\n   4: bigWig (.bw)")
    ("odir",        value<string>()->default_value("parse2wigdir"),	  "output directory name")
    ("rcenter", value<int>()->default_value(0), "consider length around the center of fragment ")
    ;
  options_description optsingle("For single-end read",100);
  optsingle.add_options()
    ("nomodel",   "predefine the fragment length (default: estimated by hamming distance plot)")
    ("flen",        value<int>()->default_value(150), "predefined fragment length\n(Automatically calculated in paired-end mode)")
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
    ("mp",        value<string>(),	  "Mappability file")
    ("mpthre",    value<double>()->default_value(0.3),	  "Threshold of low mappability regions")
    ;
  options_description optgc("GC bias normalization\n   (require large time and memory)",100);
  optgc.add_options()
    ("genome",     value<string>(),	  "reference genome sequence for GC content estimation")
    ("flen4gc",    value<int>()->default_value(120),  "fragment length for calculation of GC distribution")
    ("gcdepthoff", "do not consider depth of GC contents")
    ;
  options_description optother("Others",100);
  optother.add_options()
    ("threads,p",    value<int>()->default_value(1),  "number of threads to launch")
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
  BPRINT("Number of threads: %1%\n") % values["threads"].as<int>();
  if (!values.count("pair")) {
    cout << "Single-end mode: ";
    BPRINT("fragment length will be estimated from hamming distance\n");
    if (values.count("nomodel")) BPRINT("Predefined fragment length: %1%\n") % values["flen"].as<int>();
  } else {
    cout << "Paired-end mode: ";
    BPRINT("Maximum fragment length: %1%\n") % values["maxins"].as<int>();
  }
  if (!values.count("nofilter")) {
    BPRINT("PCR bias filtering: ON\n");
    if (values["thre_pb"].as<int>()) BPRINT("PCR bias threshold: > %1%\n") % values["thre_pb"].as<int>();
  }
  BPRINT("\t%1% reads used for library complexity\n") % values["ncmp"].as<int>();
  if (values.count("bed")) BPRINT("Bed file: %1%\n") % values["bed"].as<string>();

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
    BPRINT("\tFile directory: %1%\n") % values["mp"].as<string>();
    BPRINT("\tLow mappablitiy threshold: %1%\n") % values["mpthre"].as<double>();
  }
  if (values.count("genome")) {
    printf("Correcting GC bias:\n");
    BPRINT("\tChromosome directory: %1%\n") % values["genome"].as<string>();
  }
  printf("======================================\n");
  return;
}

void print_SeqStats(const variables_map &values, ofstream &out, const SeqStats &p, const Mapfile &mapfile)
{
  /* genome data */
  out << p.name << "\t" << p.getlen()  << "\t" << p.getlenmpbl() << "\t" << p.getpmpbl() << "\t";
  /* total reads*/
  out << boost::format("%1%\t%2%\t%3%\t%4$.1f%%\t")
    % p.bothnread() % p.seq[STRAND_PLUS].nread % p.seq[STRAND_MINUS].nread
    % (p.bothnread()*100/(double)mapfile.genome.bothnread());

  /* nonredundant reads */
  printr(out, p.bothnread_nonred(), p.bothnread());
  p.seq[STRAND_PLUS].printnonred(out);
  p.seq[STRAND_MINUS].printnonred(out);
  printr(out, p.bothnread_red(), p.bothnread());
  p.seq[STRAND_PLUS].printred(out);
  p.seq[STRAND_MINUS].printred(out);

  /* reads after GCnorm */
  if(values.count("genome")) {
    printr(out, p.bothnread_afterGC(), p.bothnread());
    p.seq[STRAND_PLUS].printafterGC(out);
    p.seq[STRAND_MINUS].printafterGC(out);
  }
  out << boost::format("%1$.3f\t") % p.depth;
  if(p.w) out << boost::format("%1$.3f\t") % p.w; else out << " - \t";
  if(values["ntype"].as<string>() == "NONE") out << p.bothnread_nonred() << "\t"; else out << p.bothnread_rpm() << "\t";

  p.gcov.print(out, mapfile.gv);
  
  p.ws.printPoispar(out);
  if(values.count("bed")) out << boost::format("%1$.3f\t") % p.FRiP;

  p.ws.printZINBpar(out);
  
  out << endl;
  return;
}

void output_stats(const variables_map &values, const Mapfile &p)
{
  string filename = p.oprefix + "." + IntToString(values["binsize"].as<int>()) + ".csv";
  ofstream out(filename);

  out << "parse2wig version " << VERSION << endl;
  out << "Input file: \"" << values["input"].as<string>() << "\"" << endl;
  out << "Redundancy threshold: >" << p.thre4filtering << endl;

  if(p.tv) out << boost::format("Library complexity: (%1$.3f) (%2%/%3%)\n") % p.complexity() % p.nt_nonred % p.nt_all;
  else     out << boost::format("Library complexity: %1$.3f (%2%/%3%)\n")   % p.complexity() % p.nt_nonred % p.nt_all;
 
  if(!values.count("nomodel")) out << "Estimated fragment length: " << p.dist.eflen << endl;
  else out << "Predefined fragment length: " << p.flen_def << endl;
  if(values.count("genome")) out << "GC summit: " << p.maxGC << endl;

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
  out << endl;
  out << "\t\t\t\t";
  out << "both\tforward\treverse\t% genome\t";
  out << "both\tforward\treverse\t";
  out << "both\tforward\treverse\t";
  if(values.count("genome")) out << "both\tforward\treverse\t";
  out << endl;

  // SeqStats
  print_SeqStats(values, out, p.genome, p);
  for(auto x:p.chr) print_SeqStats(values, out, x, p);
  
  cout << "stats is output in " << filename << "." << endl;

  return;
}

vector<char> makeGcovArray(const variables_map &values, SeqStats &chr, Mapfile &p, double r4cmp)
{
  vector<char> array;
  if(values.count("mp")) array = readMpbl_binary(values["mp"].as<string>(), ("chr" + p.lchr->name), chr.getlen());
  else array = readMpbl_binary(chr.getlen());
  if(values.count("bed")) arraySetBed(array, chr.name, p.vbed);

  int val(0);
  int size = array.size();
  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto &x: chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      
      if(rand() >= r4cmp) val=COVREAD_ALL; else val=COVREAD_NORM;
      
      int s(max(0, min(x.F3, x.F5)));
      int e(min(max(x.F3, x.F5), size-1));
      if(s >= size || e < 0) {
	cerr << "Warning: " << chr.name << " read " << s <<"-"<< e << " > array size " << array.size() << endl;
      }
      for(int i=s; i<=e; ++i) if(array[i]==MAPPABLE) array[i]=val;
    }
  }
  return array;
}

void calcGcovchr(const variables_map &values, Mapfile &p, int s, int e, double r4cmp, boost::mutex &mtx)
{
  for(int i=s; i<=e; ++i) {
    cout << p.chr[i].name << ".." << flush;
    auto array = makeGcovArray(values, p.chr[i], p, r4cmp);
    p.chr[i].gcov.calcGcov(array);
    p.genome.gcov.addGcov(p.chr[i].gcov, mtx);
  }
}

void calcGenomeCoverage(const variables_map &values, Mapfile &p)
{
  cout << "calculate genome coverage.." << flush;

  // ignore peak region
  double r = NUM_GCOV/(double)(p.genome.bothnread_nonred() - p.genome.nread_inbed);
  if(r>1){
    cerr << "Warning: number of reads is < "<< (int)(NUM_GCOV/NUM_1M) << " million.\n";
    p.gv = 1;
  }
  double r4cmp = r*RAND_MAX;

  boost::thread_group agroup;
  boost::mutex mtx;
  for(uint i=0; i<p.vsepchr.size(); i++) {
    agroup.create_thread(bind(calcGcovchr, boost::cref(values), boost::ref(p), p.vsepchr[i].s, p.vsepchr[i].e, r4cmp, boost::ref(mtx)));
  }
  agroup.join_all();
  
  cout << "done." << endl;
  return;
}

void calcFRiP(SeqStats &chr, const vector<bed> vbed)
{
  vector<char> array(chr.getlen(), MAPPABLE);
  arraySetBed(array, chr.name, vbed);
  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto &x: chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      int s(min(x.F3, x.F5));
      int e(max(x.F3, x.F5));
      for(int i=s; i<=e; ++i) {
	if(array[i]==INBED) {
	  x.inpeak = 1;
	  chr.nread_inbed++;
	  break;
	}
      }
    }
  }
  chr.setFRiP();
  return;
}

void output_wigstats(const variables_map &values, Mapfile &p)
{
  string filename = p.oprefix + "." + IntToString(values["binsize"].as<int>()) + ".binarray_dist.csv";
  ofstream out(filename);

  cout << "generate " << filename << ".." << flush;

  out << "\tGenome\t\t\t";
  for (auto x:p.chr) out << x.name << "\t\t\t\t";
  out << endl;
  out << "read number\tnum of bins genome\tprop\tZINB estimated\t";
  for (auto x:p.chr) out << "num of bins\tprop\tPoisson estimated\tZINB estimated\t";
  out << endl;

  for(size_t i=0; i<p.genome.ws.wigDist.size(); ++i) {
    out << i << "\t";
    p.genome.ws.printwigDist(out, i);
    out << p.genome.ws.getZINB(i) << "\t";
    //    out << p.genome.getZIP(i) << "\t";
    for (auto x:p.chr) {
      x.ws.printwigDist(out, i);
      out << x.ws.getPoisson(i) << "\t";
      out << x.ws.getZINB(i) << "\t";
    }
    out << endl;
  }

  cout << "done." << endl;
  return;
}
