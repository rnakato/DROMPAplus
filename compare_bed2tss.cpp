#include <string>
#include <iostream>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include "readdata.h"
#include "gene_bed.h"
#include "macro.h"

using Variables = boost::program_options::variables_map;

Variables argv_init(int argc, char* argv[])
{
  boost::program_options::options_description allopts("Options");
  allopts.add_options()
    ("genefile,g", boost::program_options::value<std::string>(), "Gene file")
    ("bed,b",      boost::program_options::value<std::string>(), "Bed file")
    ("gt",         boost::program_options::value<std::string>(), "Genome table (tab-delimited file describing the name and length of each chromosome)")
    ("mode,m",     boost::program_options::value<int>()->default_value(0),    "0: with TSS; 1: with whole gene; 2:genome coverage")
    ("updist,u",   boost::program_options::value<int>()->default_value(5000), "Allowed upstream distance from TSS")
    ("downdist,d", boost::program_options::value<int>()->default_value(5000), "Allowed downstream distance from TSS")
    ("limconv,l",  boost::program_options::value<int>()->default_value(10000), "Maxmum distance between genes for convergent sites")
    ("all,a", "Output also non-neighboring sites")
    ("nottss", "Output non-neighboring sites only")
    ("name,n", "Output name instead of id")
    ("conv,c", "Consider convergent sites")
    ("redundant", "redundant mode (output multiple genes for each peak)")
    ("refFlat,r", "refFlat format as input (default: gtf)")
    ("gene",      "Gene-level comparison (default: transcript-level)")
    ("intron,i",  "consider exon and intron (default: whole-genic)")
    ("macsxls",   "macs xls file as input")
    ("bed12", "bed12 format as input")
    ("help,h", "print this message")
    ;
  
  Variables values;
  
  if (argc==1) {
    std::cout << "\n" << allopts << std::endl;
    exit(0);
  }
  try {
    boost::program_options::parsed_options parsed = parse_command_line(argc, argv, allopts);
    store(parsed, values);
    
    if (values.count("help")) {
      std::cout << "\n" << allopts << std::endl;
      exit(0);
    }
    std::vector<std::string> opts = {};
    for (auto x: {"genefile", "bed", "gt"}) {
      if (!values.count(x)) {
	std::cerr << "specify --" << x << " option." << std::endl;
	exit(0);
      }
    }
    notify(values);

  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    exit(0);
  }
  return values;
}

template <class T>
void merge_tss2bed(const Variables &values, const HashOfGeneDataMap &mp, std::vector<T> &vbed)
{
  int updist = -values["updist"].as<int>();
  int downdist = values["downdist"].as<int>();
  for (T &x: vbed) {
    int d, dmin(-1);
    if (mp.find(x.bed.chr) != mp.end()) {
      for (auto itr = mp.at(x.bed.chr).begin(); itr != mp.at(x.bed.chr).end(); ++itr) {
	if (itr->second.strand == "+") d = x.bed.summit - itr->second.txStart;
	else                           d = itr->second.txEnd - x.bed.summit;

	if((x.gene.st != TSS && !x.gene.d) || abs(d) < dmin) {
	  dmin = abs(d);
	  x.gene.d = d;
	  x.gene.gene = &itr->second;
	}
	if (d > updist && d < downdist) {
	  x.gene.st = TSS;
	  
	  hitGene g;
	  g.st = TSS;
	  g.d = d;
	  g.gene = &itr->second;
	  x.genelist.push_back(g);
	}
      }
    }
  }

  tssdist d;
  for (auto x: vbed) d.inc(x.gene.d, x.gene.st);
  d.print();

  BPRINT("# Sites from upstream %1% bp to downstream %2% bp around TSS\n") % updist % downdist;

  vbed[0].printHead();
  std::cout << "\tfrom TSS\ttranscript name\tgene name\tstrand\ttxStart\ttxEnd" << std::endl;
  for (auto x: vbed) {
    if(values.count("all") ||
       (!values.count("nottss") && x.gene.st == TSS) ||
       (values.count("nottss")  && x.gene.st != TSS))
      x.printWithTss(values.count("redundant"));
  }

  return;
}

template <class T>
void merge_gene2bed(const Variables &values, const HashOfGeneDataMap &mp, std::vector<T> &vbed)
{
  int updist = values["updist"].as<int>();
  int downdist = values["downdist"].as<int>();

  std::vector<convsite> vconv;
  if(values.count("conv")) vconv = gen_convergent(values["limconv"].as<int>(), mp);

  for (T &x: vbed) {
    if(mp.find(x.bed.chr) == mp.end()) continue;

    int on=0;
    if(values.count("conv")) on = scanBedConv(x, vconv);
    if(on) continue;

    scanBedGene(x, mp, updist, downdist);
  }

  gdist d;
  for (auto x: vbed) d.inc(x.gene.st);
  BPRINT("# Input sites total: %1%, upstream: %2%, downstream: %3%, genic: %4%, intergenic: %5%") % d.genome % d.up % d.down % d.genic % d.inter;  
  if(values.count("conv")) BPRINT(", convergent: %1%, divergent: %2%, parallel: %3%\n") % d.conv % d.div % d.par;
  else std::cout << std::endl;
  
  vbed[0].printHead();
  std::cout << "\tstatus\ttranscript name\tgene name\tstrand\ttxStart\ttxEnd" << std::endl;
  for (auto x: vbed) x.printWithGene(values.count("redundant"));
  
  return;
}

void stUpdate(int &r, status st)
{
  if(r < st) r = st;
  return;
}

std::vector<int> makeGenomeArray(const Variables &values, std::string chr, int size,
			    const HashOfGeneDataMap &mp, std::vector<convsite> &vconv)
{
  int i,j;
  int updist = values["updist"].as<int>();
  int downdist = values["downdist"].as<int>();
  
  std::vector<int> array(size, 0);
  if(mp.find(chr) == mp.end()) goto final;

  if(values.count("conv")) {
    for (auto conv: vconv) {
      if(conv.chr == chr) {
	if( conv.end >= size) std::cout << chr <<"," << conv.chr<<"," << conv.start<<"," << conv.end <<"," << size<< std::endl;
	for(i=conv.start; i<=conv.end; i++) if(array[i] < conv.st) stUpdate(array[i], conv.st); //conv.st
      }
    }
  }

  for(auto itr = mp.at(chr).begin(); itr != mp.at(chr).end(); ++itr) {
    std::string strand = itr->second.strand;
    int s = itr->second.txStart;
    int e = itr->second.txEnd;
    //    std::cout <<itr->first << ","<< s << "-" << e << "," << strand << "," << size<< std::endl;
    
    if(strand == "+") {
      int l=std::max(0, s - updist);
      int r=std::min(e + downdist, size-1);
      for(i=l; i<=s; i++) stUpdate(array[i], UPSTREAM);
      for(i=e; i<=r; i++) stUpdate(array[i], DOWNSTREAM);
    } else {
      int l=std::max(0, s - updist);
      int r=std::min(e + updist, size-1);
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

void print_gdist(const Variables &values, gdist n, std::string str)
{
  std::cout << str << "\t";
  if(values.count("intron")) BPRINT("%1%\t%2%\t%3%\t%4%\t%5%\t%6%") % n.genome % n.up % n.down % n.exon % n.intron % n.inter;
  else BPRINT("%1%\t%2%\t%3%\t%4%\t%5%") % n.genome % n.up % n.genic % n.intron % n.inter;
  if(values.count("conv")) BPRINT("\t%1%\t%2%\t%3%\n") % n.conv % n.div % n.par;
  else std::cout << std::endl;

  std::cout << str << " (%)\t";
  if(values.count("intron")) printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 100.0, n.ratio(UPSTREAM), n.ratio(DOWNSTREAM), n.ratio(EXON), n.ratio(INTRON), n.ratio(INTERGENIC));
  else printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 100.0, n.ratio(UPSTREAM), n.ratio(DOWNSTREAM), n.ratio(GENIC), n.ratio(INTERGENIC));
  if(values.count("conv")) printf("\t%.2f\t%.2f\t%.2f\n", n.ratio(CONVERGENT), n.ratio(DIVERGENT), n.ratio(PARALLEL));
  else std::cout << std::endl;
  return;
}

template <class T>
void count_genome(const Variables &values, const HashOfGeneDataMap &mp, std::vector<T> &vbed){

  gdist n,s;
  std::vector<convsite> vconv;
  if(values.count("conv")) vconv = gen_convergent(values["limconv"].as<int>(), mp);

  auto gt = read_genometable(values["gt"].as<std::string>());
  for(auto chr:gt) {
    auto array = makeGenomeArray(values, chr.getname(), chr.getlen(), mp, vconv);
    
    for (T &x: vbed) {
      if(x.bed.chr == chr.getname()) s.inc(array[x.bed.summit]);
    }
    
    for(size_t i=0; i<array.size(); i++) n.inc(array[i]);
  }
  
  if(values.count("intron")) std::cout << "\tGenome\tupstream\tdownstream\texon\tintron\tintergenic";
  else std::cout << "\tGenome\tupstream\tdownstream\tgenic\tintergenic";
  if(values.count("conv")) std::cout << "\tconvergent\tdivergent\tparallel" << std::endl;
  else std::cout << std::endl;

  print_gdist(values, s, "peak");
  print_gdist(values, n, "base");
  
  std::cout << "relative ratio\t";
  if(values.count("intron")) printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 1.0, s.ratio(UPSTREAM)/n.ratio(UPSTREAM), s.ratio(DOWNSTREAM)/n.ratio(DOWNSTREAM), s.ratio(EXON)/n.ratio(EXON), s.ratio(INTRON)/n.ratio(INTRON), s.ratio(INTERGENIC)/n.ratio(INTERGENIC));
  else printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 100.0, s.ratio(UPSTREAM)/n.ratio(UPSTREAM), s.ratio(DOWNSTREAM)/n.ratio(DOWNSTREAM), s.ratio(GENIC)/n.ratio(GENIC), s.ratio(INTERGENIC)/n.ratio(INTERGENIC));
  if(values.count("conv")) printf("\t%.2f\t%.2f\t%.2f\n", s.ratio(CONVERGENT)/n.ratio(CONVERGENT), s.ratio(DIVERGENT)/n.ratio(DIVERGENT), s.ratio(PARALLEL)/n.ratio(PARALLEL));
  else std::cout << std::endl;
  return;
}

template <class T>
void func(const Variables &values, HashOfGeneDataMap &mp, std::vector<T> &vbed)
{
  int mode = values["mode"].as<int>();
  if(!mode) merge_tss2bed(values, mp, vbed);
  else if (mode==1) merge_gene2bed(values, mp, vbed);
  else if (mode==2) count_genome(values, mp, vbed);
  return;
}

template <class T>
void compare_bed(const Variables &values, std::string filename)
{
  auto vbed = parseBed<bed_gene<T>>(filename);
  //  printBed(vbed);

  HashOfGeneDataMap tmp;
  if(values.count("refFlat")) tmp = parseRefFlat(values["genefile"].as<std::string>());
  else                        tmp = parseGtf(values["genefile"].as<std::string>(), values.count("name"));
  
  //printMap(tmp);

  if(values.count("gene")) {
    auto gmp = construct_gmp(tmp);
    func(values, gmp, vbed);
  } else func(values, tmp, vbed);

  return;
}

int main(int argc, char* argv[])
{
  Variables values = argv_init(argc, argv);

  std::string filename(values["bed"].as<std::string>());
  if(values.count("bed12"))        compare_bed<bed12>(values, filename);
  else if(values.count("macsxls")) compare_bed<macsxls>(values, filename);
  else                             compare_bed<bed>(values, filename);

  return 0;
}

