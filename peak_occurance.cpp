#include <string>
#include <iostream>
#include <algorithm>
#include <random>
#include <boost/format.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include "readdata.h"
#include "gene_bed.h"
#include "cmdline.h"
#include "alglib.h"

using namespace boost::accumulators;
cmdline::parser argv_init(int argc, char* argv[])
{
  cmdline::parser p;

  p.add<string>("genefile", 'g', "gene file", true, "");
  p.add<string>("bed", 'b', "bed file", true, "");
  p.add<string>("genelist", 0, "gene list to be compared", true, "");
  p.add<int>("mode", 'm', "0: with TSS; 1: with whole gene", false, 0);
  p.add<int>("updist", 'u', "allowed upstream distance from TSS", false, 5000);
  p.add<int>("downdist", 'd', "allowed downstream distance from TSS", false, 5000);
  p.add<int>("limconv",  'l', "limit distance between genes for convergent sites", false, 10000);
  p.add<int>("permutation", 'p', "random permutation time", false, 1000);
  p.add("refFlat", 'r', "refFlat format as input (default: gtf)");
  p.add("macsxls", 0, "macs xls file as input");
  p.add("bed12",   0, "bed12 format as input");
  p.add("help",  'h', "print this message");
  p.add("name",  'n', "output name instead of id");
  p.add("conv",   0, "consider convergent sites");
  p.add("gene",    0, "gene-level comparison (default: transcript-level)");

  if (argc==1 || !p.parse(argc, argv) || p.exist("help")) {
    if (argc==1 || p.exist("help")) cout << "compare bed file with gene annotation file" << endl;
    cout << p.error_full() << p.usage();
    exit(1);
  }
  return p;
}

template <class T>
tssdist func(cmdline::parser p, unordered_map<string, unordered_map<string, genedata>> mp, vector<T> vbed){
  for (T &x: vbed) {
    int updist = -p.get<int>("updist");
    int downdist = p.get<int>("downdist");

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
  
  return d;
}

template <class T>
unordered_map<string, unordered_map<string, genedata>>
					    generate_rand(cmdline::parser p,
							  unordered_map<string, unordered_map<string, genedata>> mp,
							  vector<string> vgname, vector<T> &vbed, int size_emp)
{
  std::random_device random_device;  // 乱数生成器
  std::mt19937 random_engine{random_device()}; // メルセンヌ・ツイスター 32bit
  std::uniform_int_distribution<int> die_distribution{0, (int)vgname.size()-1};

  vector<string> glist_rand;
  int i(0);
  while(i<size_emp) {
    int rand = die_distribution(random_engine);
    string name = vgname[rand];
    int on(0);
    for(auto x: glist_rand) {
      if(x == name){
	on++;
	break;
      }
    }
    if(!on){
      glist_rand.push_back(name);
      i++;
    }
  }
  //  for(auto x: glist_rand) cout << x << endl;

  auto rmp = extract_mp(mp, glist_rand);
  //  countmp(rmp);
  //int n_rmp = countmp(rmp);
  //  cout << "gene num: " << n_rmp << endl;
  return rmp;
}

template <class T>
void compare_tss(cmdline::parser p,
		 unordered_map<string, unordered_map<string, genedata>> mp,
		 vector<T> &vbed, vector<string> glist)
{
  // genelist
  auto emp = extract_mp(mp, glist);
  int n_emp = countmp(emp);
  tssdist d_list = func(p, emp, vbed);

  // random
  accumulator_set<double, stats<tag::mean, tag::variance>> vn1, vn5, vn10, vn100, vnover100;
  vector<string> vgname = scanGeneName(mp);
  int max = p.get<int>("permutation");
  cout << "random permutation: " << endl;
  for(int i=0; i<max; ++i) {
    cout << i << ".." << flush;
    auto rmp = generate_rand(p, mp, vgname, vbed, n_emp);
    tssdist d_rand = func(p, rmp, vbed);
    vn1(d_rand.n1);
    vn5(d_rand.n5);
    vn10(d_rand.n10);
    vn100(d_rand.n100);
    vnover100(d_rand.nover100);
    //    d_rand.print();
  }
  //  d_list.print();
  cout << "gene num: " << n_emp << endl;  
  cout << "\t<1k\t1k~5k\t5k~10k\t10k~100k\t100k~"  << endl;
  cout << "list\t" << d_list.n1 << "\t" << d_list.n5 << "\t" << d_list.n10 << "\t" << d_list.n100 << "\t" << d_list.nover100 << endl;
  cout << "mean\t" << mean(vn1) << "\t" << mean(vn5) << "\t" << mean(vn10) << "\t" << mean(vn100) << "\t" << mean(vnover100) << endl;
  cout << "SD\t" << sqrt(variance(vn1)) << "\t" << sqrt(variance(vn5)) << "\t" << sqrt(variance(vn10)) << "\t" << sqrt(variance(vn100)) << "\t" << sqrt(variance(vnover100)) << endl;

  cout << "pvalue\t"
       << stdNormdist(d_list.n1,       mean(vn1),       sqrt(variance(vn1))) << "\t"
       << stdNormdist(d_list.n5,       mean(vn5),       sqrt(variance(vn5))) << "\t"
       << stdNormdist(d_list.n10,      mean(vn10),      sqrt(variance(vn10))) << "\t"
       << stdNormdist(d_list.n100,     mean(vn100),     sqrt(variance(vn100))) << "\t"
       << stdNormdist(d_list.nover100, mean(vnover100), sqrt(variance(vnover100))) << endl;
  return;
}

template <class T>
gdist func_gene(cmdline::parser p, unordered_map<string, unordered_map<string, genedata>> mp, vector<T> vbed)
{
  int updist = p.get<int>("updist");
  int downdist = p.get<int>("downdist");
  
  vector<convsite> vconv;
  if(p.exist("conv")) vconv = gen_convergent(p.get<int>("limconv"), mp);
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
  
  return d;
}

template <class T>
void merge_gene2bed(const cmdline::parser p,
		    const unordered_map<string, unordered_map<string, genedata>> &mp,
		    vector<T> &vbed, vector<string> glist)
{
  // genelist
  auto emp = extract_mp(mp, glist);
  int n_emp = countmp(emp);
  gdist d_list = func_gene(p, emp, vbed);

  // random
  accumulator_set<double, stats<tag::mean, tag::variance>> vup, vdown, vgenic, vinter, vconv, vdiv, vpar;
  vector<string> vgname = scanGeneName(mp);
  int max = p.get<int>("permutation");
  cout << "random permutation: " << endl;
  for(int i=0; i<max; ++i) {
    cout << i << ".." << flush;
    auto rmp = generate_rand(p, mp, vgname, vbed, n_emp);
    gdist d_rand = func_gene(p, rmp, vbed);
    vup(d_rand.up);
    vdown(d_rand.down);
    vgenic(d_rand.genic);
    vinter(d_rand.inter);
  }
  cout << "done." << endl;
  
  cout << "gene num: " << n_emp << ", peak num: " << d_list.genome << endl;
  cout << "\tupstream\tdownstream\tgenic\tintergenic";
  if(p.exist("conv")) cout << "\tconvergent\tdivergent\tparallel";
  cout << endl;
  cout << "list\t" << d_list.up << "\t" << d_list.down << "\t" << d_list.genic << "\t" << d_list.inter;
  if(p.exist("conv")) cout << "\t" << d_list.conv << "\t" << d_list.div << "\t" << d_list.par;
  cout << endl;
  cout << "mean\t" << mean(vup) << "\t" << mean(vdown) << "\t" << mean(vgenic) << "\t" << mean(vinter);
  if(p.exist("conv")) cout << "\t" << mean(vconv) << "\t" << mean(vdiv) << "\t" << mean(vpar);
  cout << endl;
  cout << "SD\t" << sqrt(variance(vup)) << "\t" << sqrt(variance(vdown)) << "\t" << sqrt(variance(vgenic)) << "\t" << sqrt(variance(vinter));
  if(p.exist("conv")) cout << "\t" << sqrt(variance(vconv)) << "\t" << sqrt(variance(vdiv)) << "\t" << sqrt(variance(vpar));
  cout << endl;
  cout << "pvalue\t"
       << stdNormdist(d_list.up,       mean(vup),      sqrt(variance(vup))) << "\t"
       << stdNormdist(d_list.down,     mean(vdown),    sqrt(variance(vdown))) << "\t"
       << stdNormdist(d_list.genic,    mean(vgenic),   sqrt(variance(vgenic))) << "\t"
       << stdNormdist(d_list.inter,    mean(vinter),   sqrt(variance(vinter)));
  if(p.exist("conv")) {
    cout << "\t" << stdNormdist(d_list.conv,   mean(vconv),    sqrt(variance(vconv)))
	 << "\t" << stdNormdist(d_list.div,    mean(vdiv),     sqrt(variance(vdiv)))
	 << "\t" << stdNormdist(d_list.par,    mean(vpar),     sqrt(variance(vpar)));
  }
  cout << endl;

  return;
}

template <class T>
void compare_bed(cmdline::parser p, string filename)
{
  auto vbed = parseBed<bed_gene<T>>(filename);
  //  printBed(vbed);

  unordered_map<string, unordered_map<string, genedata>> tmp;
  if(p.exist("refFlat")) tmp = parseRefFlat(p.get<string>("genefile"));
  else                   tmp = parseGtf(p.get<string>("genefile"),  p.exist("name"));
  
  //printMap(tmp);
  
  auto glist = readGeneList(p.get<string>("genelist"));

  if(p.exist("gene")) {
    auto gmp = construct_gmp(tmp);
    if(!p.get<int>("mode")) compare_tss(p, gmp, vbed, glist);
    else if (p.get<int>("mode")==1) merge_gene2bed(p, gmp, vbed, glist);
  } else {
    if(!p.get<int>("mode")) compare_tss(p, tmp, vbed, glist);
    else if (p.get<int>("mode")==1) merge_gene2bed(p, tmp, vbed, glist);
  }
  
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

