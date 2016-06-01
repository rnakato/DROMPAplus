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
#include <seqan/bam_io.h>
#include <seqan/store.h>

using namespace boost::program_options;
using namespace seqan;

class FragmentSingle {
public:
  CharString name;
  int chr_F3;
  int F3;
  Strand strand;
  int fraglen;
  int readlen_F3;
  //  short num_multimapped;
  FragmentSingle(const BamAlignmentRecord &record, int flen=0):
    name(record.qName),
    chr_F3(record.rID),
    fraglen(flen),
    readlen_F3(length(record.seq))
  {
    if(!fraglen) fraglen = abs(record.tLen);
    if(hasFlagRC(record)) {
      strand = STRAND_MINUS;
      F3 = record.beginPos + length(record.seq) -1;
    } else {
      strand = STRAND_PLUS;
      F3 = record.beginPos;
    }
  }
  void print() {
    cout << name << ", " << chr_F3 << ", " << F3<< ", "<< strand << ", " << "fraglen " << fraglen << "," <<readlen_F3 << endl;
  }
};

class FragmentPair: public FragmentSingle {
public:
  int chr_F5;
  int F5;
  int readlen_F5;
  //  short num_multimapped;
  FragmentPair(const BamAlignmentRecord &record):
    FragmentSingle(record)
  {}
};

void printRecord(BamAlignmentRecord record)
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
}

vector<FragmentSingle> do_bampe(const variables_map &values, BamFileIn & bamFileIn)
{
  vector<FragmentSingle> vFrag;
  int nreads(0);
  int maxins(values["maxins"].as<int>());
  BamAlignmentRecord record;
  while(!atEnd(bamFileIn)) {
    readRecord(record, bamFileIn);
    //      printRecord(record);
    if(!hasFlagAllProper(record) || hasFlagLast(record)) continue;
    FragmentSingle frag(record);
    //    frag.print();
    if(frag.fraglen > maxins) continue;
    nreads++;
    vFrag.push_back(frag);
  }
  cout << "loaded " << nreads << " mapped reads\n" << endl;
  return vFrag;
}

vector<FragmentSingle> do_bamse(const variables_map &values, BamFileIn & bamFileIn)
{
  vector<FragmentSingle> vFrag;
  int nreads(0);
  int flen = values["flen"].as<int>();
  BamAlignmentRecord record;
  while(!atEnd(bamFileIn)) {
    readRecord(record, bamFileIn);
    //      printRecord(record);
    
    if(hasFlagUnmapped(record) || hasFlagQCNoPass(record)) continue;
    if(hasFlagFirst(record) || hasFlagLast(record)) cerr << "Warning: parsing paired-end file as single-end." << endl;
    FragmentSingle frag(record, flen);  ////// fraglenを推定することにすれば？？？？
    //    frag.print();
    nreads++;
    vFrag.push_back(frag);
  }
  cout << "loaded " << nreads << " mapped reads\n" << endl;
  return vFrag;
}

void parse_sam(const variables_map &values, string inputfile, RefGenome &g)
{
  CharString filename(inputfile);
  BamFileIn bamFileIn;
  if (!open(bamFileIn, toCString(filename))) {
    cerr << "ERROR: Could not open " << filename << endl;
    exit(1);
  }

  try {
    // useFragStore(bamFileIn);
    BamHeader header;
    readHeader(header, bamFileIn);
    
    vector<FragmentSingle> vFrag;
    if (values.count("pair")) vFrag = do_bampe(values, bamFileIn);
    else vFrag = do_bamse(values, bamFileIn);
    cout << "loaded " << vFrag.size() << " mapped reads\n" << endl;
  }
  catch (Exception const & e) {
    cout << "ERROR: " << e.what() << endl;
    exit(1);
  }
  return;
}

void do_parse(const variables_map &values, string inputfile, RefGenome &g){
  isFile(inputfile);
  BPRINT("Parsing %1%...\n") % inputfile;

  string ftype = values["ftype"].as<string>();
  if(ftype == "SAM" || ftype == "BAM") parse_sam(values, inputfile, g);
  //else if(ftype == "BOWTIE")   parse_bowtie(values, inputfile, g); 
  // else if(ftype == "TAGALIGN") parse_tagAlign(values, inputfile, g);

  printf("done.\n");
}

void read_mapfile(const variables_map &values, RefGenome &g){
  vector<string> v;
  boost::split(v, values["input"].as<string>(), boost::algorithm::is_any_of(","));
  for(auto x:v) do_parse(values, x, g);
  
  /*  Mapfile *mapfile = mapfile_new(g->chrnum, p);
  
  add_SeqStats_to_genome(mapfile, g);
  output_read_fragment_distribution(p, mapfile);*/  /* output distributions of read length and fragment length */

  return;
}
