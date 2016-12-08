#include <string>
#include <iostream>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include "SSP/src/readdata.h"
#include "gene_bed.h"

using namespace std;
using namespace boost::program_options;

variables_map argv_init(int argc, char* argv[])
{
  options_description allopts("Options");
  allopts.add_options()
    ("genefile,g", value<string>(), "Gene file")
    ("bed,b",      value<vector<string>>(), "Bed file")
    ("gt",         value<string>(), "Genome table (tab-delimited file describing the name and length of each chromosome)")
    ("updist,u",   value<int>()->default_value(5000), "Allowed upstream distance from TSS")
    ("downdist,d", value<int>()->default_value(5000), "Allowed downstream distance from TSS")
    ("name,n", "Output name instead of id")
    ("refFlat,r", "refFlat format as input (default: gtf)")
    ("length,l",  "output summed length (default: number of peaks)")
    ("gene",      "Gene-level comparison (default: transcript-level)")
    ("macsxls",   "macs xls file as input")
    ("bed12", "bed12 format as input")
    ("help,h", "print this message")
    ;
  
  variables_map values;
  
  if (argc==1) {
    cout << "\n" << allopts << endl;
    exit(0);
  }
  try {
    parsed_options parsed = parse_command_line(argc, argv, allopts);
    store(parsed, values);
    
    if (values.count("help")) {
      cout << "\n" << allopts << endl;
      exit(0);
    }
    vector<string> opts = {};
    for (auto x: {"genefile", "bed"}) {
      if (!values.count(x)) {
	cerr << "specify --" << x << " option." << endl;
	exit(0);
      }
    }
    notify(values);

  } catch (exception &e) {
    cout << e.what() << endl;
    exit(0);
  }
  return values;
}

template <class T>
void merge_tss2bed(const variables_map &values,
		   const unordered_map<string, unordered_map<string, genedata>> &mp,
		   vector<vector<T>> &vbedlist)
{
  int updist = -values["updist"].as<int>();
  int downdist = values["downdist"].as<int>();
  
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
	    if(values.count("length")) nbed[i] += x.length();
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
void compare_bed(const variables_map &values)
{
  vector<vector <T>> vbedlist;
  for(auto x: values["bed"].as<vector<string>>()) {
    //    cout << x << endl;
    auto vbed = parseBed<T>(x);
    //    printBed(vbed);
    vbedlist.push_back(vbed);
  }

  auto merge_bed(vbedlist);

  unordered_map<string, unordered_map<string, genedata>> tmp; // hash for transcripts
  if(values.count("refFlat")) tmp = parseRefFlat(values["genefile"].as<string>());
  else                        tmp = parseGtf(values["genefile"].as<string>());

  //printMap(tmp);

  if(values.count("gene")) {
    auto gmp = construct_gmp(tmp);              // hash for genes
    merge_tss2bed(values, gmp, vbedlist);
  } else merge_tss2bed(values, tmp, vbedlist);
  
  return;
}

int main(int argc, char* argv[])
{
  variables_map values = argv_init(argc, argv);

  if(values.count("bed12"))        compare_bed<bed12>(values);
  else if(values.count("macsxls")) compare_bed<macsxls>(values);
  else                             compare_bed<bed>(values);

  return 0;
}

