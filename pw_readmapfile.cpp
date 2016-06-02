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
#include <seqan/store.h>
#include <ext/stdio_filebuf.h>

using namespace boost::program_options;
using namespace seqan;

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

/*void do_bampe(const variables_map &values, Mapfile &p, BamFileIn & bamFileIn)
{
  int maxins(values["maxins"].as<int>());
  BamAlignmentRecord record;
  while(!atEnd(bamFileIn)) {
    readRecord(record, bamFileIn);
    //      printRecord(record);
    if(!hasFlagAllProper(record) || hasFlagLast(record)) continue;
    FragmentSingle frag(record);
    //frag.print();
    if(frag.fraglen > maxins) continue;
    p.addfrag(frag);
  }
  return;
  }*/

template <class T>
void do_bamse(const variables_map &values, Mapfile &p, T & in)
{
  int flen = values["flen"].as<int>();
 
  string lineStr; 
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0]=='@') continue;

    cout << lineStr << endl;
    exit(0);
    vector<string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    FragmentSingle frag(v);////// fraglenを推定することにすれば？？？？

  }

  /*
  BamAlignmentRecord record;
  while(!atEnd(bamFileIn)) {
    readRecord(record, bamFileIn);
    //      printRecord(record);
    if(hasFlagUnmapped(record) || hasFlagQCNoPass(record)) continue;
    if(hasFlagFirst(record) || hasFlagLast(record))
      cerr << "Warning: parsing paired-end file as single-end." << endl;
    FragmentSingle frag(record, flen);////// fraglenを推定することにすれば？？？？
    //    frag.print();
    p.addfrag(frag);
    }*/
  return;
}

void parse_sam(const variables_map &values, string inputfile, Mapfile &p, RefGenome &g)
{
  //  CharString filename(inputfile);
  //BamFileIn bamFileIn;
  //if (!open(bamFileIn, toCString(filename))) PRINTERR("Could not open " << filename);

  try {

    if(values["ftype"].as<string>()=="SAM") {
      ifstream in(inputfile);
      if(!in) PRINTERR("Could not open " << inputfile << ".");
      do_bamse<ifstream>(values, p, in);      
//      if (values.count("pair")) do_bampe(values, p, in);
      // else do_bamse(values, p, in);
    } else if(values["ftype"].as<string>()=="BAM") {
      string command = "samtools view -h " + inputfile;
      FILE *fp = popen(command.c_str(), "r");
      __gnu_cxx::stdio_filebuf<char> *p_fb = new __gnu_cxx::stdio_filebuf<char>(fp, ios_base::in);
      istream in(static_cast<streambuf *>(p_fb));
      do_bamse<istream>(values, p, in);


      /*      char command[100];
      char str[256];
      sprintf(command, "samtools view -h %s", inputfile.c_str());
      FILE *IN = popen(command, "r");
      while((fgets(str, 256, IN))!=NULL){ 
	if(str[0]=='\n' || str[0]=='@') continue;
	printf("%s\n",str);
	}*/
      //      if (values.count("pair")) do_bampe(values, p, in);
      //else do_bamse(values, p, in);
    }
    // useFragStore(bamFileIn);
    //    BamHeader header;
    //readHeader(header, bamFileIn);

  }
  catch (Exception const & e) {
    PRINTERR(e.what());
  }
  return;
}

void read_mapfile(const variables_map &values, RefGenome &g){
  Mapfile p(g);
  //  vector<FragmentSingle> vFrag;

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
