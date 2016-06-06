/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <algorithm>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <ext/stdio_filebuf.h>
#include "pw_gv.h"
#include "pw_readmapfile.h"
#include "util.h"
#include "common.h"
#include "macro.h"

using namespace boost::program_options;

void check_redundant_reads(const variables_map &values, Mapfile &p);
void filtering_eachchr_single(const variables_map &values, Mapfile &p, strandData &seq);
void filtering_eachchr_pair(const variables_map &values, Mapfile &p, SeqStats &chr);

template <class T>
void do_bampe(const variables_map &values, Mapfile &p, T &in)
{
  int maxins(values["maxins"].as<int>());

  string lineStr; 
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0]=='@') continue;

    vector<string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    int sv = stoi(v[1]); // bitwise FLAG
    if(sv&4 || sv&512 || sv&1024) continue;
    if(!(sv&2)) continue;
    if(sv&128) {
      p.addF5(v[9].length());
      continue;
    }
    FragmentSingle frag(v);
    if(frag.fraglen > maxins) continue;
    //frag.print();
    p.addfrag(frag);
  }

  return;
}

template <class T>
void do_bamse(const variables_map &values, Mapfile &p, T & in)
{
  string lineStr; 
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0]=='@') continue;
    vector<string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    int sv = stoi(v[1]); // bitwise FLAG
    // unmapped reads, low quality reads
    if(sv&4 || sv&512 || sv&1024) continue;
    if(sv&64 || sv&128) cerr << "Warning: parsing paired-end file as single-end." << endl;
    FragmentSingle frag(v, p.dist.eflen);  ////// fraglenを推定することにすれば？？？？
    //    frag.print();
    p.addfrag(frag);
  }

  return;
}

void parse_sam(const variables_map &values, string inputfile, Mapfile &p)
{
  if(values["ftype"].as<string>()=="SAM") {  // SAM
    ifstream in(inputfile);
    if(!in) PRINTERR("Could not open " << inputfile << ".");
    if (values.count("pair")) do_bampe(values, p, in);
    else do_bamse(values, p, in);
  }
  else if(values["ftype"].as<string>()=="BAM") {  // BAM
    string command = "samtools view -h " + inputfile;
    FILE *fp = popen(command.c_str(), "r");
    __gnu_cxx::stdio_filebuf<char> *p_fb = new __gnu_cxx::stdio_filebuf<char>(fp, ios_base::in);
    istream in(static_cast<streambuf *>(p_fb));
    if (values.count("pair")) do_bampe(values, p, in);
    else do_bamse(values, p, in);
  }

  return;
}

void outputDist(const variables_map &values, Mapfile &p){
  string outputfile = values["odir"].as<string>() + "/" + values["output"].as<string>() + ".readlength_dist.xls";

  p.dist.printReadlen(outputfile, values.count("pair"));
  if(values.count("pair")) {
    outputfile = values["odir"].as<string>() + "/" + values["output"].as<string>() + ".fragmentlength_dist.xls";
    p.dist.printFraglen(outputfile);
  }

  return;
}

void read_mapfile(const variables_map &values, Mapfile &p){
  vector<string> v;
  boost::split(v, values["input"].as<string>(), boost::algorithm::is_any_of(","));
  for(auto inputfile: v) {
    isFile(inputfile);
    BPRINT("Parsing %1%...\n") % inputfile;
    string ftype = values["ftype"].as<string>();
    if(ftype == "SAM" || ftype == "BAM") parse_sam(values, inputfile, p);
    //else if(ftype == "BOWTIE")   parse_bowtie(values, inputfile, g); 
    // else if(ftype == "TAGALIGN") parse_tagAlign(values, inputfile, g);
    printf("done.\n");
  }

  /* output distributions of read length and fragment length */
  outputDist(values, p);

  /* PCR bias filtering and ignore enrichregions */
  check_redundant_reads(values, p);
  p.update();

  // calculate depth
  p.calcdepth();

#ifdef DEBUG
  p.printstats();
#endif

  return;
}

// for single-end
void hashFilterAll(unordered_map<int, int> &mp, strandData &seq, const int thre)
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

void hashFilterCmp(unordered_map<int, int> &mp, Mapfile &p, const strandData &seq, const int thre)
{
  for(auto x: seq.vRead){
    if(rand() >= p.r4cmp) continue;
    p.nt_all++;
    if(mp.find(x.F3) != mp.end()) {
      if(mp[x.F3] < thre) {
	++mp[x.F3];
	p.nt_nonred++;
      } else {
	p.nt_red++;
      }
    } else {
      mp[x.F3] = 1;
      p.nt_nonred++;
    }
  }
  return;
}

// for paired-end
void hashFilterAll(unordered_map<string, int> &mp, strandData &seq, const int thre)
{
  for(auto &x:seq.vRead) {
    int Fmin = min(x.F3, x.F5);
    int Fmax = max(x.F3, x.F5);
    string str = IntToString(Fmin) + "-" + IntToString(Fmax);
    if(mp.find(str) != mp.end()) {
      if(mp[str] < thre) {
	++mp[str];
	seq.nread_nonred++;
      } else {
	x.duplicate = 1;
	seq.nread_red++;
      }
    } else {
      mp[str] = 1;
      seq.nread_nonred++;
    }
  }
  return;
}

void hashFilterCmp(unordered_map<string, int> &mp, Mapfile &p, const strandData &seq, const int thre)
{
  for(auto x: seq.vRead){
    if(rand() >= p.r4cmp) continue;
    int Fmin = min(x.F3, x.F5);
    int Fmax = max(x.F3, x.F5);
    string str = IntToString(Fmin) + "-" + IntToString(Fmax);
    p.nt_all++;
    if(mp.find(str) != mp.end()) {
      if(mp[str] < thre) {
	++mp[str];
	p.nt_nonred++;
      } else {
	p.nt_red++;
      }
    } else {
      mp[str] = 1;
      p.nt_nonred++;
    }
  }
  return;
}

void check_redundant_reads(const variables_map &values, Mapfile &p)
{
  int threshold;
  int strand;
  if(values["thre_pb"].as<int>()) threshold = values["thre_pb"].as<int>();
  else threshold = max(1, (int)(p.genome.bothnread() *10/(double)p.genome.len_mpbl));
  
  cout << "\nChecking redundant reads: redundancy threshold " << threshold << endl;
  p.thre4filtering = threshold;

  // Library complexity
  double r = values["ncmp"].as<int>()/(double)p.genome.bothnread();
  if(r>1){
    cerr << "Warning: number of reads is < "<< (int)(values["ncmp"].as<int>()/NUM_1M) <<" million.\n";
    p.tv = 1;
  }
  p.r4cmp = r*RAND_MAX;

  for (auto &chr: p.chr) {
     cout << chr.name << ".." << flush;
     if (values.count("pair")) filtering_eachchr_pair(values, p, chr);
     else {
       for(strand=0; strand<STRANDNUM; strand++) filtering_eachchr_single(values, p, chr.seq[strand]);
     }
  }

#ifdef DEBUG
  BPRINT("nt_all %1%, nt_nonred %2%, nt_red %3%, complexity %4%\n") % p.nt_all % p.nt_nonred % p.nt_red % p.complexity();
#endif

  printf("done.\n");
  return;
}


void filtering_eachchr_single(const variables_map &values, Mapfile &p, strandData &seq)
{
  unordered_map<int, int> mp;
  hashFilterAll(mp, seq, p.thre4filtering);
  
#ifdef DEBUG
  BPRINT("nread: %1%, nonred: %2%, red: %3%\n") % seq.nread % seq.nread_nonred % seq.nread_red;
#endif

  unordered_map<int, int> mp2;
  hashFilterCmp(mp2, p, seq, p.thre4filtering);

  return;
}

void filtering_eachchr_pair(const variables_map &values, Mapfile &p, SeqStats &chr)
{
  //  int Fmin, Fmax;
  unordered_map<string, int> mp;
  for(int strand=0; strand<STRANDNUM; ++strand) {
    hashFilterAll(mp, chr.seq[strand], p.thre4filtering);
  }

#ifdef DEBUG
  for(int strand=0; strand<STRANDNUM; strand++) 
    BPRINT("%4% strand %5% nread: %1%, nonred: %2%, red: %3%\n") % chr.seq[strand].nread % chr.seq[strand].nread_nonred % chr.seq[strand].nread_red % chr.name % strand;
#endif

  unordered_map<string, int> mp2;
  for(int strand=0; strand<STRANDNUM; strand++){
    hashFilterCmp(mp2, p, chr.seq[strand], p.thre4filtering);
  }

  return;
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


//#define SEQAN_HAS_ZLIB 1
//#include <zlib.h>
//#include <seqan/store.h>
//using namespace seqan;
/*void printRecord(BamAlignmentRecord record)
{
  cout << "AllProper " << hasFlagAllProper(record) << endl;
  cout << "Duplicate " << hasFlagDuplicate(record) << endl;
  cout << "First " << hasFlagFirst(record)     << endl;
  cout << "Last " << hasFlagLast(record)      << endl;
  cout << "Multiple " << hasFlagMultiple(record)  << endl;
  cout << "NextRC " << hasFlagNextRC(record)    << endl;
  cout << "NextUnmapped " << hasFlagNextUnmapped(record) << endl;
  cout << "QCNoPass " << hasFlagQCNoPass(record)  << endl;
  cout << "RC " << hasFlagRC(record)        << endl;
  cout << "Secondary " << hasFlagSecondary(record) << endl;
  cout << "Unmapped " << hasFlagUnmapped(record)  << endl;
}

void useFragStore(BamFileIn bamFileIn){
  FragmentStore<> fragStore;
  readRecords(fragStore, bamFileIn);
  uint nreads = length(fragStore.alignedReadStore);
  cout << "loaded " << nreads << " mapped reads\n" << endl;

  for (uint i=0; i<nreads; ++i) {
    cout << fragStore.alignedReadStore[i].beginPos<< endl;
    cout << fragStore.alignedReadStore[i].contigId<< endl;
    cout << fragStore.alignedReadStore[i].endPos<< endl;
    cout << fragStore.alignedReadStore[i].id<< endl;
    //      cout << fragStore.alignedReadStore[i].pairMatchId << endl;
    //cout << fragStore.alignedReadStore[i].INVALID_ID << endl;
    cout << fragStore.alignedReadStore[i].readId<< endl;
    //      FragmentData flag(fragStore.alignedReadStore[i]);
  }
  return;
  }*/


  //  CharString filename(inputfile);
  //BamFileIn bamFileIn;
  //if (!open(bamFileIn, toCString(filename))) PRINTERR("Could not open " << filename);

  // useFragStore(bamFileIn);
  //    BamHeader header;
  //readHeader(header, bamFileIn);
  
