#include <string>
#include <iostream>
#include <algorithm>
#include <boost/format.hpp>
#include "readdata.h"
#include "gene_bed.h"
#include "warn.h"
#include "cmdline.h"

cmdline::parser argv_init(int argc, char* argv[])
{
  cmdline::parser p;
  p.add<string>("genefile", 'g', "gene file", true, "");
  p.add<string>("bed", 'b', "bed file", true, "");
  p.add<string>("gt",  0, "genome table", false, "");
  p.add<int>("mode", 'm', "0: with TSS; 1: with whole gene; 2:genome coverage", false, 0);
  p.add<int>("updist", 'u', "allowed upstream distance from TSS", false, 5000);
  p.add<int>("downdist", 'd', "allowed downstream distance from TSS", false, 5000);
  p.add<int>("limconv",  'l', "limit distance between genes for convergent sites", false, 10000);
  p.add("all",   'a', "output also non-neighboring sites");
  p.add("name",  'n', "output name instead of id");
  p.add("conv",   'c', "consider convergent sites");
  p.add("refFlat", 'r', "refFlat format as input (default: gtf)");
  p.add("gene",    0, "gene-level comparison (default: transcript-level)");
  p.add("intron",  'i', "consider exon and intron (default: whole-genic)");
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
		    vector<T> &vbed)
{
  int updist = -p.get<int>("updist");
  int downdist = p.get<int>("downdist");
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
    if((x.st == TSS || p.exist("all"))) x.printWithTss();
  }

  return;
}

template <class T>
void merge_gene2bed(const cmdline::parser p,
		    const unordered_map<string, unordered_map<string, genedata>> &mp,
		    vector<T> &vbed)
{
  int updist = p.get<int>("updist");
  int downdist = p.get<int>("downdist");

  vector<convsite> vconv;
  if(p.exist("conv")) vconv = gen_convergent(p, mp);
  int convflag = p.exist("conv");
  
  for (T &x: vbed) {
    if(mp.find(x.bed.chr) == mp.end()) continue;

    int on=0;
    if(convflag) on = scanBedConv(x, vconv);
    if(on) continue;
    
    scanBedGene(x, mp, updist, downdist);
  }

  gdist d;
  for (auto x: vbed) d.inc(x.st);
  cout << boost::format("# Input sites total: %1%, upstream: %2%, downstream: %3%, genic: %4%, intergenic: %5%") % d.genome % d.up % d.down % d.genic % d.inter;  
  if(p.exist("conv")) cout << boost::format(", convergent: %1%, divergent: %2%, parallel: %3%\n") % d.conv % d.div % d.par;
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

vector<int> makeGenomeArray(const cmdline::parser p, string chr, int size,
			    const unordered_map<string, unordered_map<string, genedata>> &mp,
			    vector<convsite> &vconv)
{
  int i,j;
  int updist = p.get<int>("updist");
  int downdist = p.get<int>("downdist");
  
  vector<int> array(size, 0);
  if(mp.find(chr) == mp.end()) goto final;

  if(p.exist("conv")) {
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
    
    if(p.exist("intron")) {
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

template <class T>
void count_genome(const cmdline::parser p, const unordered_map<string, unordered_map<string, genedata>> &mp, vector<T> &vbed){

  gdist n,s;
  vector<convsite> vconv;
  if(p.exist("conv")) vconv = gen_convergent(p, mp);

  auto gt = read_genometable(p.get<string>("gt"));
  for(auto itr = gt.begin(); itr != gt.end(); ++itr){
    string chr(itr->first);

    //    cout << itr->first << "\t" << itr->second << endl;
    
    auto array = makeGenomeArray(p, chr, itr->second, mp, vconv);
    
    for (T &x: vbed) {
      //      if(chr == "chr19")cout << x.bed.chr <<"," << chr << endl;
      if(x.bed.chr != chr) continue;
      //     cout << x.bed.chr <<","<< x.bed.summit<<"," << chr <<","<<  itr->second << endl;
      s.inc(array[x.bed.summit]);
    }
    
    for(int i=0; i<itr->second; i++) n.inc(array[i]);
  }
  
  if(p.exist("intron")) cout << "\tGenome\tupstream\tdownstream\texon\tintron\tintergenic";
  else cout << "\tGenome\tupstream\tdownstream\tgenic\tintergenic";
  if(p.exist("conv")) cout << "\tconvergent\tdivergent\tparallel" << endl;
  else cout << endl;

  print_gdist(p, s, "peak");
  print_gdist(p, n, "base");
  
  cout << "relative ratio\t";
  if(p.exist("intron")) printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 1.0, s.ratio(UPSTREAM)/n.ratio(UPSTREAM), s.ratio(DOWNSTREAM)/n.ratio(DOWNSTREAM), s.ratio(EXON)/n.ratio(EXON), s.ratio(INTRON)/n.ratio(INTRON), s.ratio(INTERGENIC)/n.ratio(INTERGENIC));
  else printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 100.0, s.ratio(UPSTREAM)/n.ratio(UPSTREAM), s.ratio(DOWNSTREAM)/n.ratio(DOWNSTREAM), s.ratio(GENIC)/n.ratio(GENIC), s.ratio(INTERGENIC)/n.ratio(INTERGENIC));
  if(p.exist("conv")) printf("\t%.2f\t%.2f\t%.2f\n", s.ratio(CONVERGENT)/n.ratio(CONVERGENT), s.ratio(DIVERGENT)/n.ratio(DIVERGENT), s.ratio(PARALLEL)/n.ratio(PARALLEL));
  else cout << endl;
  return;
}

template <class T>
void func(cmdline::parser p, unordered_map<string, unordered_map<string, genedata>> &mp, vector<T> &vbed)
{
  if(!p.get<int>("mode")) merge_tss2bed(p, mp, vbed);
  else if (p.get<int>("mode")==1) merge_gene2bed(p, mp, vbed);
  else if (p.get<int>("mode")==2) count_genome(p, mp, vbed);
  return;
}

template <class T>
void compare_bed(cmdline::parser p, string filename)
{
  auto vbed = parseBed<bed_gene<T>>(filename);
  //  printBed(vbed);

  unordered_map<string, unordered_map<string, genedata>> tmp; // hash for transcripts
  if(p.exist("refFlat")) tmp = parseRefFlat(p.get<string>("genefile"));
  else                   tmp = parseGtf(p.get<string>("genefile"), p.exist("name"));
  
  //printMap(tmp);

  if(p.exist("gene")) {
    auto gmp = construct_gmp(tmp);              // hash for genes
    func(p, gmp, vbed);
  } else func(p, tmp, vbed);

  return;
}

int main(int argc, char* argv[])
{
  cmdline::parser p = argv_init(argc, argv);

  string filename = p.get<string>("bed");
  if(p.exist("bed12"))        compare_bed<bed12>(p, filename);
  else if(p.exist("macsxls")) compare_bed<macsxls>(p, filename);
  else                        compare_bed<bed>(p, filename);

  return 0;
}

