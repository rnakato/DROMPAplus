#include <string>
#include <iostream>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include "readdata.h"
#include "gene_bed.h"

using namespace boost::program_options;

variables_map argv_init(int argc, char* argv[])
{
  options_description allopts("Options");
  allopts.add_options()
    ("genefile,g", value<string>(), "Gene file")
    ("bed,b",      value<string>(), "Bed file")
    ("gt",         value<string>(), "Genome table (tab-delimited file describing the name and length of each chromosome)")
    ("mode,m",     value<int>()->default_value(0),    "0: with TSS; 1: with whole gene; 2:genome coverage")
    ("updist,u",   value<int>()->default_value(5000), "Allowed upstream distance from TSS")
    ("downdist,d", value<int>()->default_value(5000), "Allowed downstream distance from TSS")
    ("limconv,l",  value<int>()->default_value(10000), "Maxmum distance between genes for convergent sites")
    ("all,a", "Output also non-neighboring sites")
    ("name,n", "Output name instead of id")
    ("conv,c", "Consider convergent sites")
    ("refFlat,r", "refFlat format as input (default: gtf)")
    ("gene",      "Gene-level comparison (default: transcript-level)")
    ("intron,i",  "consider exon and intron (default: whole-genic)")
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
    for (auto x: {"genefile", "bed", "gt"}) {
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
		    vector<T> &vbed)
{
  int updist = -values["updist"].as<int>();
  int downdist = values["downdist"].as<int>();
  for (T &x: vbed) {
    int d, dmin(-1);
    if(mp.find(x.bed.chr) != mp.end()) {
      for(auto itr = mp.at(x.bed.chr).begin(); itr != mp.at(x.bed.chr).end(); ++itr) {
	if(itr->second.strand == "+") d = x.bed.summit - itr->second.txStart;
	else                          d = itr->second.txEnd - x.bed.summit;

	if((x.st != TSS && !x.d) || dmin > abs(d)){
	  dmin = abs(d);
	  x.d = d;
	  x.gene = &itr->second;
	}
	if(d > updist && d < downdist) x.st = TSS;
      }
    }
  }

  tssdist d;
  for (auto x: vbed) d.inc(x.d, x.st);
  d.print();

  cout << boost::format("# Sites from upstream %1% bp to downstream %2% bp around TSS\n") % updist % downdist;

  vbed[0].printHead();
  cout << "\tfrom TSS\ttranscript name\tgene name\tstrand\ttxStart\ttxEnd" << endl;
  for (auto x: vbed) {
    if(x.st == TSS || values.count("all")) x.printWithTss();
  }

  return;
}

template <class T>
void merge_gene2bed(const variables_map &values,
		    const unordered_map<string, unordered_map<string, genedata>> &mp,
		    vector<T> &vbed)
{
  int updist = values["updist"].as<int>();
  int downdist = values["downdist"].as<int>();

  vector<convsite> vconv;
  if(values.count("conv")) vconv = gen_convergent(values["limconv"].as<int>(), mp);
  
  for (T &x: vbed) {
    if(mp.find(x.bed.chr) == mp.end()) continue;

    int on=0;
    if(values.count("conv")) on = scanBedConv(x, vconv);
    if(on) continue;
    
    scanBedGene(x, mp, updist, downdist);
  }

  gdist d;
  for (auto x: vbed) d.inc(x.st);
  cout << boost::format("# Input sites total: %1%, upstream: %2%, downstream: %3%, genic: %4%, intergenic: %5%") % d.genome % d.up % d.down % d.genic % d.inter;  
  if(values.count("conv")) cout << boost::format(", convergent: %1%, divergent: %2%, parallel: %3%\n") % d.conv % d.div % d.par;
  else cout << endl;
  
  vbed[0].printHead();
  cout << "\tstatus\ttranscript name\tgene name\tstrand\ttxStart\ttxEnd" << endl;
  for (auto x: vbed) x.printWithGene();
  
  return;
}

void stUpdate(int &r, status st)
{
  if(r < st) r = st;
  return;
}

vector<int> makeGenomeArray(const variables_map &values, string chr, int size,
			    const unordered_map<string, unordered_map<string, genedata>> &mp,
			    vector<convsite> &vconv)
{
  int i,j;
  int updist = values["updist"].as<int>();
  int downdist = values["downdist"].as<int>();
  
  vector<int> array(size, 0);
  if(mp.find(chr) == mp.end()) goto final;

  if(values.count("conv")) {
    for (auto conv: vconv) {
      if(conv.chr == chr) {
	if( conv.end >= size) cout << chr <<"," << conv.chr<<"," << conv.start<<"," << conv.end <<"," << size<< endl;
	for(i=conv.start; i<=conv.end; i++) if(array[i] < conv.st) stUpdate(array[i], conv.st); //conv.st
      }
    }
  }

  for(auto itr = mp.at(chr).begin(); itr != mp.at(chr).end(); ++itr) {
    string strand = itr->second.strand;
    int s = itr->second.txStart;
    int e = itr->second.txEnd;
    //    cout <<itr->first << ","<< s << "-" << e << "," << strand << "," << size<< endl;
    
    if(strand == "+") {
      int l=max(0, s - updist);
      int r=min(e + downdist, size-1);
      for(i=l; i<=s; i++) stUpdate(array[i], UPSTREAM);
      for(i=e; i<=r; i++) stUpdate(array[i], DOWNSTREAM);
    } else {
      int l=max(0, s - updist);
      int r=min(e + updist, size-1);
      for(i=e; i<=r; i++) stUpdate(array[i], UPSTREAM);
      for(i=l; i<=s; i++) stUpdate(array[i], DOWNSTREAM);
    }
    
    if(values.count("intron")) {
      for(i=0; i<itr->second.exonCount; i++) {
	int es = itr->second.exon[i].start;
	int ee = itr->second.exon[i].end;
	for(j=es; j<=ee; j++) stUpdate(array[j], EXON);
	if(i < itr->second.exonCount-1) {
	  es = itr->second.exon[i].end +1;
	  ee = itr->second.exon[i+1].start;
	  for(j=es; j<ee; j++) stUpdate(array[j], INTRON);
	}
      }
    } else {
      for(i=s; i<=e; i++) stUpdate(array[i], GENIC);
    }
  }
  
 final:
  return array;
}

void print_gdist(const variables_map &values, gdist n, string str)
{
  cout << str << "\t";
  if(values.count("intron")) cout << boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%") % n.genome % n.up % n.down % n.exon % n.intron % n.inter;
  else cout << boost::format("%1%\t%2%\t%3%\t%4%\t%5%") % n.genome % n.up % n.genic % n.intron % n.inter;
  if(values.count("conv")) cout << boost::format("\t%1%\t%2%\t%3%\n") % n.conv % n.div % n.par;
  else cout << endl;

  cout << str << " (%)\t";
  if(values.count("intron")) printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 100.0, n.ratio(UPSTREAM), n.ratio(DOWNSTREAM), n.ratio(EXON), n.ratio(INTRON), n.ratio(INTERGENIC));
  else printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 100.0, n.ratio(UPSTREAM), n.ratio(DOWNSTREAM), n.ratio(GENIC), n.ratio(INTERGENIC));
  if(values.count("conv")) printf("\t%.2f\t%.2f\t%.2f\n", n.ratio(CONVERGENT), n.ratio(DIVERGENT), n.ratio(PARALLEL));
  else cout << endl;
  return;
}

template <class T>
void count_genome(const variables_map &values, const unordered_map<string, unordered_map<string, genedata>> &mp, vector<T> &vbed){

  gdist n,s;
  vector<convsite> vconv;
  if(values.count("conv")) vconv = gen_convergent(values["limconv"].as<int>(), mp);

  auto gt = read_genometable(values["gt"].as<string>());
  for(auto itr = gt.begin(); itr != gt.end(); ++itr){
    string chr(itr->first);

    //    cout << itr->first << "\t" << itr->second << endl;
    
    auto array = makeGenomeArray(values, chr, itr->second, mp, vconv);
    
    for (T &x: vbed) {
      if(x.bed.chr == chr) s.inc(array[x.bed.summit]);
    }
    
    for(int i=0; i<itr->second; i++) n.inc(array[i]);
  }
  
  if(values.count("intron")) cout << "\tGenome\tupstream\tdownstream\texon\tintron\tintergenic";
  else cout << "\tGenome\tupstream\tdownstream\tgenic\tintergenic";
  if(values.count("conv")) cout << "\tconvergent\tdivergent\tparallel" << endl;
  else cout << endl;

  print_gdist(values, s, "peak");
  print_gdist(values, n, "base");
  
  cout << "relative ratio\t";
  if(values.count("intron")) printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 1.0, s.ratio(UPSTREAM)/n.ratio(UPSTREAM), s.ratio(DOWNSTREAM)/n.ratio(DOWNSTREAM), s.ratio(EXON)/n.ratio(EXON), s.ratio(INTRON)/n.ratio(INTRON), s.ratio(INTERGENIC)/n.ratio(INTERGENIC));
  else printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 100.0, s.ratio(UPSTREAM)/n.ratio(UPSTREAM), s.ratio(DOWNSTREAM)/n.ratio(DOWNSTREAM), s.ratio(GENIC)/n.ratio(GENIC), s.ratio(INTERGENIC)/n.ratio(INTERGENIC));
  if(values.count("conv")) printf("\t%.2f\t%.2f\t%.2f\n", s.ratio(CONVERGENT)/n.ratio(CONVERGENT), s.ratio(DIVERGENT)/n.ratio(DIVERGENT), s.ratio(PARALLEL)/n.ratio(PARALLEL));
  else cout << endl;
  return;
}

template <class T>
void func(const variables_map &values, unordered_map<string, unordered_map<string, genedata>> &mp, vector<T> &vbed)
{
  int mode = values["mode"].as<int>();
  if(!mode) merge_tss2bed(values, mp, vbed);
  else if (mode==1) merge_gene2bed(values, mp, vbed);
  else if (mode==2) count_genome(values, mp, vbed);
  return;
}

template <class T>
void compare_bed(const variables_map &values, string filename)
{
  auto vbed = parseBed<bed_gene<T>>(filename);
  //  printBed(vbed);

  unordered_map<string, unordered_map<string, genedata>> tmp; // hash for transcripts
  if(values.count("refFlat")) tmp = parseRefFlat(values["genefile"].as<string>());
  else                        tmp = parseGtf(values["genefile"].as<string>(), values.count("name"));
  
  //printMap(tmp);

  if(values.count("gene")) {
    auto gmp = construct_gmp(tmp);              // hash for genes
    func(values, gmp, vbed);
  } else func(values, tmp, vbed);

  return;
}

int main(int argc, char* argv[])
{
  variables_map values = argv_init(argc, argv);

  string filename(values["bed"].as<string>());
  if(values.count("bed12"))        compare_bed<bed12>(values, filename);
  else if(values.count("macsxls")) compare_bed<macsxls>(values, filename);
  else                             compare_bed<bed>(values, filename);

  return 0;
}

