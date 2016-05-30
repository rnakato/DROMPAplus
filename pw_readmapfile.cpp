/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <boost/algorithm/string.hpp>
#include <seqan/bam_io.h>
#include <seqan/store.h>
#include "pw_gv.h"
#include "pw_readmapfile.h"
#include "util.h"
#include "macro.h"

using namespace boost::program_options;
using namespace seqan;

/*class FragmentData{
public:
  string name;
  int chr_F3, chr_F5;
  int F3, F5;
  Strand strand;
  int fraglen;
  int readlen_F3, readlen_F5;
  //  short num_multimapped;
  FragmentData(const BamAlignmentRecord &record):
    name(record.qName), chr_F3(record.rID)
    //    readlen_F3(length(record.seq))
    //    fraglen(record.tLen)
  {
    if(hasFlagRC(record)) {
      strand = STRAND_MINUS;
      F3 = record.beginPos + length(record.seq) -1;
    } else {
      strand = STRAND_PLUS;
      F3 = record.beginPos;
    }
  }
  };*/

void parse_sam(const variables_map &values, string inputfile, RefGenome &g)
{
  CharString filename(inputfile);
  BamFileIn bamFileIn;
  if (!open(bamFileIn, toCString(filename))) {
    cerr << "ERROR: Could not open " << filename << endl;
    exit(1);
  }

  FragmentStore<> fragStore;
  //  int nreads(0);
  try {
    // Copy header
    /* BamHeader header;
       readHeader(header, bamFileIn);
       int chr_num=length(header.sequenceInfos);
       for(int i=0; i<chr_num; ++i) {
       cout << header.sequenceInfos[i].i1 << endl; // name
       cout << header.sequenceInfos[i].i2 << endl; // length
     }*/
    readRecords(fragStore, bamFileIn);
    uint nreads = length(fragStore.alignedReadStore);
    cout << nreads << endl;
    
    /*    BamAlignmentRecord record;
    while(!atEnd(bamFileIn)) {
      readRecord(record, bamFileIn);
      
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
      
      if(hasFlagUnmapped(record) || hasFlagQCNoPass(record)) continue;
      FragmentData flag(record);
      if(frag.fraglen > values["maxins"].as<int>() && frag.fraglen < 0) continue;
      nreads++;
      }*/
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
