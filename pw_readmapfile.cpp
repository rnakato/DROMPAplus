/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "pw_readmapfile.h"
#include <algorithm>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <ext/stdio_filebuf.h>
#include <time.h>
#include "pw_shiftprofile.h"

void parseSam(const MyOpt::Variables &, const std::string &, Mapfile &);
void parseBowtie(const MyOpt::Variables &, const std::string &, Mapfile &);
void parseTagAlign(const MyOpt::Variables &, const std::string &, Mapfile &);
void filtering_eachchr_single(const MyOpt::Variables &values, Mapfile &, SeqStats &);
void filtering_eachchr_pair(const MyOpt::Variables &values, Mapfile &, SeqStats &);

void read_mapfile(const MyOpt::Variables &values, Mapfile &p)
{
  std::vector<std::string> v;
  boost::split(v, values["input"].as<std::string>(), boost::algorithm::is_any_of(","));
  for(auto inputfile: v) {
    isFile(inputfile);
    BPRINT("Parsing %1%...\n") % inputfile;
    std::string ftype = values["ftype"].as<std::string>();
    if(ftype == "SAM" || ftype == "BAM") parseSam(values, inputfile, p);
    else if(ftype == "BOWTIE") parseBowtie(values, inputfile, p);
    else if(ftype == "TAGALIGN") parseTagAlign(values, inputfile, p);
    printf("done.\n");
  }
  p.genome.setnread();

  if(!p.genome.bothnread()) PRINTERR("no read in input file.");

  p.setFraglen(values);
  
  /* output distributions of read length and fragment length */
  p.outputDistFile(values);

  return;
}

template <class T>
void do_bampe(const MyOpt::Variables &values, Mapfile &p, T &in)
{
  int maxins(values["maxins"].as<int>());

  std::string lineStr;
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0]=='@') continue;

    std::vector<std::string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    int sv(stoi(v[1]));   // bitwise FLAG
    if(sv&4 || sv&512 || sv&1024) continue;
    if(!(sv&2)) continue;
    if(sv&128) {
      p.addF5(v[9].length());
      continue;
    }
    Fragment frag;
    frag.addSAM(v, values.count("pair"), sv);
    if(frag.fraglen > maxins) continue;
    //frag.print();
    p.addfrag(frag);
  }

  return;
}

template <class T>
void do_bamse(const MyOpt::Variables &values, Mapfile &p, T & in)
{
  std::string lineStr; 
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0]=='@') continue;
    std::vector<std::string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    int sv(stoi(v[1])); // bitwise FLAG
    // unmapped reads, low quality reads
    if(sv&4 || sv&512 || sv&1024) continue;
    if(sv&64 || sv&128) std::cerr << "Warning: parsing paired-end file as single-end." << std::endl;
    Fragment frag;
    frag.addSAM(v, values.count("pair"), sv);
    //    std::cout << lineStr << std::endl;
    // frag.print();
    p.addfrag(frag);
  }

  return;
}

void parseSam(const MyOpt::Variables &values, const std::string &inputfile, Mapfile &p)
{
  if(values["ftype"].as<std::string>()=="SAM") {  // SAM
    std::ifstream in(inputfile);
    if(!in) PRINTERR("Could not open " << inputfile << ".");
    if (values.count("pair")) do_bampe(values, p, in);
    else do_bamse(values, p, in);
  }
  else if(values["ftype"].as<std::string>()=="BAM") {  // BAM
    std::string command = "samtools view -h " + inputfile;
    FILE *fp = popen(command.c_str(), "r");
    __gnu_cxx::stdio_filebuf<char> *p_fb = new __gnu_cxx::stdio_filebuf<char>(fp, std::ios_base::in);
    std::istream in(static_cast<std::streambuf *>(p_fb));
    if (values.count("pair")) do_bampe(values, p, in);
    else do_bamse(values, p, in);
  }

  return;
}

void parseBowtie(const MyOpt::Variables &values, const std::string &inputfile, Mapfile &p)
{
  int maxins(values["maxins"].as<int>());
  std::ifstream in(inputfile);
  if(!in) PRINTERR("Could not open " << inputfile << ".");

  std::string chr_F3(""), chr_F5(""), nametemp("");
  int F5(0);
  Fragment fragpair;
  
  std::string lineStr;
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty()) continue;

    std::vector<std::string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));

    if (values.count("pair")) {
      std::vector<std::string> read;
      boost::split(read, v[0], boost::algorithm::is_any_of("/"));
      if(nametemp != "" && nametemp != read[0]) PRINTERR("Error:  Invalid read pair." << nametemp <<"-" << read[0]);
      if(read[1] == "1") {  // F3 read
	chr_F3 = rmchr(v[2]);
	fragpair.readlen_F3 = v[4].length();
	if(v[1] == "+") { 
	  fragpair.strand = STRAND_PLUS;
	  fragpair.F3 = stoi(v[3]);
	} else {
	  fragpair.strand = STRAND_MINUS;
	  fragpair.F3 = stoi(v[3]) + fragpair.readlen_F3;
	}
      } else {  
	chr_F5 = rmchr(v[2]);
	if(v[1] == "+") F5 = stoi(v[3]);
	else            F5 = stoi(v[3]) + v[4].length();
	p.addF5(v[4].length());
      }
      if(chr_F3 != "" && chr_F5 != ""){
	if(chr_F3 == chr_F5) {
	  fragpair.chr = chr_F3;
	  fragpair.fraglen = abs(F5 - fragpair.F3);
	  if(fragpair.fraglen <= maxins) p.addfrag(fragpair);
	  //	  fragpair.print();
	}
	chr_F3 = "";
	chr_F5 = "";
	nametemp = "";
      }
    } else {
      if(v[0].find("/2") != std::string::npos) PRINTERR("Warning: parsing paired-end file as single-end");
      Fragment frag;
      frag.chr = rmchr(v[2]);
      frag.readlen_F3 = v[4].length();
      if(v[1] == "+") { 
	frag.strand = STRAND_PLUS;
	frag.F3 = stoi(v[3]);
      } else {
	frag.strand = STRAND_MINUS;
	frag.F3 = stoi(v[3]) + frag.readlen_F3;
      }
      //      std::cout << lineStr << std::endl;
      //      frag.print();
      p.addfrag(frag);
    }

  }
  return;
}

template <class T>
void funcTagAlign(const MyOpt::Variables &values, Mapfile &p, T &in)
{
  std::string lineStr;
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty()) continue;

    std::vector<std::string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    if(v.size() < 6) PRINTERR("Use tagAlign (BED3+3) file");
    
    if (values.count("pair")) PRINTERR("tagAlign format does not support paired-end file.\n");
    else {
      int start(stoi(v[1]));
      int end(stoi(v[2]));
      Fragment frag;
      frag.chr = rmchr(v[0]);
      frag.readlen_F3 = abs(end - start);
      if(v[5] == "+") {
	frag.strand = STRAND_PLUS;
	frag.F3 = start;
      } else {
	frag.strand = STRAND_MINUS;
	frag.F3 = start + frag.readlen_F3;
      }
      //      std::cout << lineStr << std::endl;
      //frag.print();
      p.addfrag(frag);
    }
  }
  return;
}

void parseTagAlign(const MyOpt::Variables &values, const std::string &inputfile, Mapfile &p)
{
  if(inputfile.find(".gz") != std::string::npos) {
    std::string command = "zcat " + inputfile;
    FILE *fp = popen(command.c_str(), "r");
    __gnu_cxx::stdio_filebuf<char> *p_fb = new __gnu_cxx::stdio_filebuf<char>(fp, std::ios_base::in);
    std::istream in(static_cast<std::streambuf *>(p_fb));
    funcTagAlign(values, p, in);
  } else {
    std::ifstream in(inputfile);
    if(!in) PRINTERR("Could not open " << inputfile << ".");
    funcTagAlign(values, p, in);
  }
  return;
}

void printDist(std::ofstream &out, const std::vector<int> v, const std::string str, const long nread)
{
  out << "\n" << str << " length distribution" << std::endl;
  out << "length\tnumber\tproportion" << std::endl;
  for(size_t i=0; i<v.size(); ++i)
    if(v[i]) out << boost::format("%1%\t%2%\t%3%\n") % i % v[i] % (v[i]/static_cast<double>(nread));
  return;
}

void hashFilterAllSingle(std::unordered_map<int, int> &mp, strandData &seq, const int thre)
{
  for(auto &x:seq.vRead) {
    if(mp.find(x.F3) != mp.end()) {
      if(mp[x.F3] < thre) {
	++mp[x.F3];
	seq.nread_nonred++;
      } else {
	x.duplicate = 1;
	seq.nread_red++;
      }
    } else {
      mp[x.F3] = 1;
      seq.nread_nonred++;
    }
  }
  return;
}

void hashFilterCmpSingle(std::unordered_map<int, int> &mp, Mapfile &p, const strandData &seq, const int thre)
{
  for(auto x: seq.vRead){
    if(rand() >= p.getr4cmp()) continue;
    p.incNtAll();
    if(mp.find(x.F3) != mp.end()) {
      if(mp[x.F3] < thre) {
	++mp[x.F3];
	p.incNtNonred();
      } else {
	p.incNtRed();
      }
    } else {
      mp[x.F3] = 1;
      p.incNtNonred();
    }
  }
  return;
}

void hashFilterAllPair(std::unordered_map<std::string, int> &mp, strandData &seq, const int thre)
{
  for(auto &x:seq.vRead) {
    int Fmin = std::min(x.F3, x.F5);
    int Fmax = std::max(x.F3, x.F5);
    std::string str = IntToString(Fmin) + "-" + IntToString(Fmax);
    //    std::cout << str << std::endl;
    if(mp.find(str) != mp.end()) {
      if(mp[str] < thre) {
	++mp[str];
	++seq.nread_nonred;
      } else {
	x.duplicate = 1;
	++seq.nread_red;
      }
    } else {
      mp[str] = 1;
      ++seq.nread_nonred;
    }
  }
  return;
}

void hashFilterCmpPair(std::unordered_map<std::string, int> &mp, Mapfile &p, const strandData &seq, const int thre)
{
  for(auto x: seq.vRead){
    if(rand() >= p.getr4cmp()) continue;
    int Fmin = std::min(x.F3, x.F5);
    int Fmax = std::max(x.F3, x.F5);
    std::string str = IntToString(Fmin) + "-" + IntToString(Fmax);
    p.incNtAll();
    if(mp.find(str) != mp.end()) {
      if(mp[str] < thre) {
	++mp[str];
	p.incNtNonred();
      } else {
	p.incNtRed();
      }
    } else {
      mp[str] = 1;
      p.incNtNonred();
    }
  }
  return;
}

void checkRedundantReads(const MyOpt::Variables &values, Mapfile &p)
{
  p.setthre4filtering(values);
  
  // Library complexity
  double r = values["ncmp"].as<int>()/static_cast<double>(p.genome.bothnread());
  if(r>1){
    std::cerr << "Warning: number of reads is < "<< (int)(values["ncmp"].as<int>()/NUM_1M) <<" million.\n";
    p.lackOfRead4Complexity_on();
  }
  p.setr4cmp(r*RAND_MAX);
  
  for(uint i=0; i<p.genome.chr.size(); ++i) {
    if (values.count("pair")) filtering_eachchr_pair(values, p, p.genome.chr[i]);
    else                      filtering_eachchr_single(values, p, p.genome.chr[i]);
  }
  
  printf("done.\n");
  return;
}

void filtering_eachchr_single(const MyOpt::Variables &values, Mapfile &p, SeqStats &chr)
{
  for(int strand=0; strand<STRANDNUM; strand++) {

    std::unordered_map<int, int> mp;
    hashFilterAllSingle(mp, chr.seq[strand], p.getthre4filtering());
    
    std::unordered_map<int, int> mp2;
    hashFilterCmpSingle(mp2, p, chr.seq[strand], p.getthre4filtering());
  }
  
  return;
}

void filtering_eachchr_pair(const MyOpt::Variables &values, Mapfile &p, SeqStats &chr)
{
  std::unordered_map<std::string, int> mp;
  for(int strand=0; strand<STRANDNUM; ++strand) {
    hashFilterAllPair(mp, chr.seq[strand], p.getthre4filtering());
  }

  std::unordered_map<std::string, int> mp2;
  for(int strand=0; strand<STRANDNUM; strand++) {
    hashFilterCmpPair(mp2, p, chr.seq[strand], p.getthre4filtering());
  }

  return;
}

void estimateFragLength(const MyOpt::Variables &values, Mapfile &p)
{
  if(!values.count("pair") && !values.count("nomodel")) {
    clock_t t1 = clock();
    //    strShiftProfile(values, p, "exjaccard");
    clock_t t2 = clock();
    std::cout << "Jaccard Vec: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";

    // thresholdが2以上の時にbitを使うと、total readがおかしくなるので
        // background uniformityが1を超える可能性がある
    strShiftProfile(values, p, "jaccard"); 
    clock_t t3 = clock();
    std::cout << "Jaccard Bit: " << static_cast<double>(t3 - t2) / CLOCKS_PER_SEC << "sec.\n";

    strShiftProfile(values, p, "fvp"); 
    clock_t t4 = clock();
    std::cout << "Fragment variability: " << static_cast<double>(t4 - t3) / CLOCKS_PER_SEC << "sec.\n";
   
    
    //strShiftProfile(values, p, "hdp");
    clock_t t5 = clock();
    std::cout << "Hamming: " << static_cast<double>(t5 - t4) / CLOCKS_PER_SEC << "sec.\n";

    //    strShiftProfile(values, p, "ccp");
    clock_t t6 = clock();
    std::cout << "ccp: " << static_cast<double>(t6 - t5) / CLOCKS_PER_SEC << "sec.\n";
  }
  
  exit(0);
  
  p.genome.calcdepthAll(p.getflen(values));
  p.genome.setF5All(p.getflen(values));
}

int check_sv(int sv)
{
  // for paired-end
  /*  LOG("   the read is paired in sequencing: %d\n",sv&1);
  LOG("   the read is mapped in a proper pair: %d\n",sv&2);
  LOG("   the query sequence itself is unmapped: %d\n",sv&4);
  LOG("   the mate is unmapped: %d\n",sv&8);
  LOG("   strand of the query (1 for reverse): %d\n",sv&16);
  LOG("   strand of the mate: %d\n",sv&32);
  LOG("   the read is the first read(F3) in a pair: %d\n",sv&64);
  LOG("   the read is the second read(F5) in a pair: %d\n",sv&128);
  LOG("   the alignment is not primary: %d\n",sv&256);
  LOG("   the read fails platform/vendor quality checks: %d\n",sv&512);
  LOG("   the read is either a PCR or an optical duplicate: %d\n",sv&1024);*/

  /*  LOG("   template having multiple segments in sequencing: %d\n",sv&1);
  LOG("   each segment properly aligned according to the aligner: %d\n",sv&2);
  LOG("   segment unmapped: %d\n",sv&4);
  LOG("   next segment in the template unmapped: %d\n",sv&8);
  LOG("   SEQ being reverse complemented: %d\n",sv&16);
  LOG("   SEQ of the next segment in the template being reversed: %d\n",sv&32);
  LOG("   the first segment in the template: %d\n",sv&64);
  LOG("   the last segment in the template: %d\n",sv&128);
  LOG("   secondary alignment: %d\n",sv&256);
  LOG("   not passing quality controls: %d\n",sv&512);
  LOG("   PCR or optical duplicate: %d\n",sv&1024);
  LOG("   supplementary alignment: %d\n",sv&2048);
  */

  // unmapped reads
  if(sv&4) goto err;
  // low quality reads
  if(sv&512 || sv&1024) goto err;
  //  if(p->rtype==READTYPE_PAIR){
    // unproper pair
    if(!(sv&2)) goto err;
    // unmatched pairs and interchromosomal pairs
    if(sv&8) goto err;
    // read pair mapped in same strand (for paired-end)
    if((sv&16 && sv&32) || (!(sv&16) && !(sv&32))) goto err;
    // }

 err:
  return 0;
}
