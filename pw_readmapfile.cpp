/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */

#define SEQAN_HAS_ZLIB 1

#include <boost/algorithm/string.hpp>
#include "pw_gv.h"
#include "pw_readmapfile.h"
#include "util.h"
#include "macro.h"
#include <zlib.h>
//#include <seqan/store.h>
#include <ext/stdio_filebuf.h>

using namespace boost::program_options;
//using namespace seqan;

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
    if(!(sv&2) || sv&128) continue;
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
  int flen = values["flen"].as<int>();
 
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

    FragmentSingle frag(v, flen);  ////// fraglenを推定することにすれば？？？？
    //    frag.print();
    p.addfrag(frag);
  }

  return;
}

void parse_sam(const variables_map &values, string inputfile, Mapfile &p, RefGenome &g)
{
  if(values["ftype"].as<string>()=="SAM") {
    ifstream in(inputfile);
    if(!in) PRINTERR("Could not open " << inputfile << ".");
    do_bamse<ifstream>(values, p, in);
    if (values.count("pair")) do_bampe(values, p, in);
    else do_bamse(values, p, in);
  } else if(values["ftype"].as<string>()=="BAM") {
    string command = "samtools view -h " + inputfile;
    FILE *fp = popen(command.c_str(), "r");
    __gnu_cxx::stdio_filebuf<char> *p_fb = new __gnu_cxx::stdio_filebuf<char>(fp, ios_base::in);
    istream in(static_cast<streambuf *>(p_fb));
    if (values.count("pair")) do_bampe(values, p, in);
    else do_bamse(values, p, in);
  }

  return;
}

void read_mapfile(const variables_map &values, RefGenome &g){
  Mapfile p(g);

  vector<string> v;
  boost::split(v, values["input"].as<string>(), boost::algorithm::is_any_of(","));
  for(auto inputfile: v) {
    isFile(inputfile);
    BPRINT("Parsing %1%...\n") % inputfile;
    string ftype = values["ftype"].as<string>();
    if(ftype == "SAM" || ftype == "BAM") parse_sam(values, inputfile, p, g);
    //else if(ftype == "BOWTIE")   parse_bowtie(values, inputfile, g); 
    // else if(ftype == "TAGALIGN") parse_tagAlign(values, inputfile, g);
    printf("done.\n");
  }

  p.update();

#ifdef DEBUG
  for (auto x:p.chr) cout<< " loaded " << x.bothnread() << " mapped reads." << endl;
  cout<< " loaded " << p.nread() << " mapped reads." << endl;
  p.dist.print();
#endif

  /*   output_read_fragment_distribution(p, mapfile);*/  /* output distributions of read length and fragment length */

  return;
}

int check_sv(int sv)
{
  // for paired-end
  LOG("   the read is paired in sequencing: %d\n",sv&1);
  LOG("   the read is mapped in a proper pair: %d\n",sv&2);
  LOG("   the query sequence itself is unmapped: %d\n",sv&4);
  LOG("   the mate is unmapped: %d\n",sv&8);
  LOG("   strand of the query (1 for reverse): %d\n",sv&16);
  LOG("   strand of the mate: %d\n",sv&32);
  LOG("   the read is the first read(F3) in a pair: %d\n",sv&64);
  LOG("   the read is the second read(F5) in a pair: %d\n",sv&128);
  LOG("   the alignment is not primary: %d\n",sv&256);
  LOG("   the read fails platform/vendor quality checks: %d\n",sv&512);
  LOG("   the read is either a PCR or an optical duplicate: %d\n",sv&1024);
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
  
