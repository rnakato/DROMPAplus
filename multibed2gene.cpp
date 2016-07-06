#include <string>
#include <iostream>
#include <algorithm>
#include <boost/format.hpp>
#include "readdata.h"
#include "gene_bed.h"
#include "cmdline.h"

cmdline::parser argv_init(int argc, char* argv[])
{
  cmdline::parser p;
  p.add<string>("genefile", 'g', "gene file", true, "");
  p.add<string>("gt",  0, "genome table", false, "");
  p.add<int>("updist", 'u', "allowed upstream distance from TSS", false, 5000);
  p.add<int>("downdist", 'd', "allowed downstream distance from TSS", false, 5000);
  p.add("name",  'n', "output name instead of id");
  p.add("refFlat", 'r', "refFlat format as input (default: gtf)");
  p.add("length",  'l', "output summed length (default: number of peaks)");
  p.add("gene",    0, "gene-level comparison (default: transcript-level)");
  p.add("macsxls", 0, "macs xls file as input");
  p.add("bed12",   0, "bed12 format as input");
  p.add("help",  'h', "print this message");

  if (argc==1 || !p.parse(argc, argv) || p.exist("help")) {
    if (argc==1 || p.exist("help")) cout << "compare bed file with gene annotation file" << endl;
    cout << p.error_full() << p.usage();
    exit(1);
  }
  return p;
}

template <class T>
void merge_tss2bed(const cmdline::parser p,
		    const unordered_map<string, unordered_map<string, genedata>> &mp,
		   vector<vector<T>> &vbedlist)
{
  int updist = -p.get<int>("updist");
  int downdist = p.get<int>("downdist");
  int plen = p.exist("length");
  
  for(auto itr = mp.begin(); itr != mp.end(); ++itr) {
    string chr = itr->first;
    for(auto itr2 = mp.at(chr).begin(); itr2 != mp.at(chr).end(); ++itr2) {
      vector<int> nbed(vbedlist.size(), 0);
      string strand = itr2->second.strand;
      int s = itr2->second.txStart;
      int e = itr2->second.txEnd;
      for(uint i=0; i<vbedlist.size(); ++i) {
	for (T &x: vbedlist[i]) {
	  int d(0);
	  if(strand == "+") d = x.summit - s;
	  else              d = e - x.summit;
	  if(d > updist && d < downdist){
	    if(plen) nbed[i] += x.length();
	    else nbed[i]++;
	  }   
	}
      }
      cout << itr2->first<< "\t";
      for(uint i=0; i<vbedlist.size(); ++i) {
	cout << "\t" << nbed[i];
      }
      cout << endl;
    }
  }

  return;
}

template <class T>
void compare_bed(cmdline::parser p)
{
  vector<string> filename;
  for(size_t i=0; i<p.rest().size(); i++) filename.push_back(p.rest()[i]);

  vector<vector <T>> vbedlist;
  for(auto x: filename) {
    //    cout << x << endl;
    auto vbed = parseBed<T>(x);
    //    printBed(vbed);
    vbedlist.push_back(vbed);
  }

  auto merge_bed(vbedlist);
  
  unordered_map<string, unordered_map<string, genedata>> tmp; // hash for transcripts
  if(p.exist("refFlat")) tmp = parseRefFlat(p.get<string>("genefile"));
  else                   tmp = parseGtf(p.get<string>("genefile"),  p.exist("name"));

  //printMap(tmp);

  if(p.exist("gene")) {
    auto gmp = construct_gmp(tmp);              // hash for genes
    merge_tss2bed(p, gmp, vbedlist);
  } else merge_tss2bed(p, tmp, vbedlist);
  
  return;
}

int main(int argc, char* argv[])
{
  cmdline::parser p = argv_init(argc, argv);

  if(p.exist("bed12"))        compare_bed<bed12>(p);
  else if(p.exist("macsxls")) compare_bed<macsxls>(p);
  else                        compare_bed<bed>(p);

  return 0;
}

