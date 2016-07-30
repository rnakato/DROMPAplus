/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <algorithm>
#include <random>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include "readdata.h"
#include "util.h"
#include "gene_bed.h"
#include "alglib.h"

using namespace std;
using namespace boost::accumulators;
using namespace boost::program_options;

variables_map argv_init(int argc, char* argv[])
{
  options_description allopts("Options");
  allopts.add_options()
    ("genefile,g", value<string>(), "Gene file")
    ("bed,b",      value<string>(), "Bed file")
    ("genelist",   value<string>(), "Gene list to be compared")
    ("mode,m",     value<int>()->default_value(0),    "0: with TSS; 1: with whole gene")
    ("updist,u",   value<int>()->default_value(5000), "Allowed upstream distance from TSS")
    ("downdist,d", value<int>()->default_value(5000), "Allowed downstream distance from TSS")
    ("limconv,l",  value<int>()->default_value(10000), "Maxmum distance between genes for convergent sites")
    ("permutation,p", value<int>()->default_value(1000), "random permutation time")
    ("refFlat,r", "refFlat format as input (default: gtf)")
    ("macsxls",   "macs xls file as input")
    ("bed12", "bed12 format as input")
    ("name,n", "Output name instead of id")
    ("conv,c", "Consider convergent sites")
    ("gene",      "Gene-level comparison (default: transcript-level)")
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
    for (auto x: {"genefile", "bed", "genelist"}) {
      if (!values.count(x)) {
	cerr << "specify --" << x << " option." << endl;
	exit(0);
      }
      isFile(values[x].as<string>());
    }
    notify(values);

  } catch (exception &e) {
    cout << e.what() << endl;
    exit(0);
  }
  return values;
}

template <class T>
tssdist func(const variables_map &values, unordered_map<string, unordered_map<string, genedata>> mp, vector<T> vbed){
  for (T &x: vbed) {
    int updist = -values["updist"].as<int>();
    int downdist = values["downdist"].as<int>();
    
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
					    generate_rand(const variables_map &values,
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
void compare_tss(const variables_map &values,
		 unordered_map<string, unordered_map<string, genedata>> mp,
		 vector<T> &vbed, vector<string> glist)
{
  // genelist
  auto emp = extract_mp(mp, glist);
  int n_emp = countmp(emp);
  tssdist d_list = func(values, emp, vbed);

  // random
  accumulator_set<double, stats<tag::mean, tag::variance>> vn1, vn5, vn10, vn100, vnover100;
  vector<string> vgname = scanGeneName(mp);
  int max = values["permutation"].as<int>();
  cerr << "random permutation: " << endl;
  for(int i=0; i<max; ++i) {
    cerr << i << ".." << flush;
    auto rmp = generate_rand(values, mp, vgname, vbed, n_emp);
    tssdist d_rand = func(values, rmp, vbed);
    vn1(d_rand.n1);
    vn5(d_rand.n5);
    vn10(d_rand.n10);
    vn100(d_rand.n100);
    vnover100(d_rand.nover100);
    //    d_rand.print();
  }
  cerr << "done." << endl;
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
gdist func_gene(const variables_map &values, unordered_map<string, unordered_map<string, genedata>> mp, vector<T> vbed)
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
  
  return d;
}

template <class T>
void merge_gene2bed(const variables_map &values,
		    const unordered_map<string, unordered_map<string, genedata>> &mp,
		    vector<T> &vbed, vector<string> glist)
{
  // genelist
  auto emp = extract_mp(mp, glist);
  int n_emp = countmp(emp);
  gdist d_list = func_gene(values, emp, vbed);

  // random
  accumulator_set<double, stats<tag::mean, tag::variance>> vup, vdown, vgenic, vinter, vconv, vdiv, vpar;
  vector<string> vgname = scanGeneName(mp);
  int max = values["permutation"].as<int>();
  cerr << "random permutation: " << endl;
  for(int i=0; i<max; ++i) {
    cerr << i << ".." << flush;
    auto rmp = generate_rand(values, mp, vgname, vbed, n_emp);
    gdist d_rand = func_gene(values, rmp, vbed);
    vup(d_rand.up);
    vdown(d_rand.down);
    vgenic(d_rand.genic);
    vinter(d_rand.inter);
  }
  cerr << "done." << endl;
  
  cout << "gene num: " << n_emp << ", peak num: " << d_list.genome << endl;
  cout << "\tupstream\tdownstream\tgenic\tintergenic";
  if(values.count("conv")) cout << "\tconvergent\tdivergent\tparallel";
  cout << endl;
  cout << "list\t" << d_list.up << "\t" << d_list.down << "\t" << d_list.genic << "\t" << d_list.inter;
  if(values.count("conv")) cout << "\t" << d_list.conv << "\t" << d_list.div << "\t" << d_list.par;
  cout << endl;
  cout << "mean\t" << mean(vup) << "\t" << mean(vdown) << "\t" << mean(vgenic) << "\t" << mean(vinter);
  if(values.count("conv")) cout << "\t" << mean(vconv) << "\t" << mean(vdiv) << "\t" << mean(vpar);
  cout << endl;
  cout << "SD\t" << sqrt(variance(vup)) << "\t" << sqrt(variance(vdown)) << "\t" << sqrt(variance(vgenic)) << "\t" << sqrt(variance(vinter));
  if(values.count("conv")) cout << "\t" << sqrt(variance(vconv)) << "\t" << sqrt(variance(vdiv)) << "\t" << sqrt(variance(vpar));
  cout << endl;
  cout << "pvalue\t"
       << stdNormdist(d_list.up,       mean(vup),      sqrt(variance(vup))) << "\t"
       << stdNormdist(d_list.down,     mean(vdown),    sqrt(variance(vdown))) << "\t"
       << stdNormdist(d_list.genic,    mean(vgenic),   sqrt(variance(vgenic))) << "\t"
       << stdNormdist(d_list.inter,    mean(vinter),   sqrt(variance(vinter)));
  if(values.count("conv")) {
    cout << "\t" << stdNormdist(d_list.conv,   mean(vconv),    sqrt(variance(vconv)))
	 << "\t" << stdNormdist(d_list.div,    mean(vdiv),     sqrt(variance(vdiv)))
	 << "\t" << stdNormdist(d_list.par,    mean(vpar),     sqrt(variance(vpar)));
  }
  cout << endl;

  return;
}

template <class T>
void compare_bed(const variables_map &values, string filename)
{
  auto vbed = parseBed<bed_gene<T>>(filename);
  //  printBed(vbed);

  unordered_map<string, unordered_map<string, genedata>> tmp;
  if(values.count("refFlat")) tmp = parseRefFlat(values["genefile"].as<string>());
  else                        tmp = parseGtf(values["genefile"].as<string>(), values.count("name"));
  
  //printMap(tmp);
  
  auto glist = readGeneList(values["genelist"].as<string>());

  int mode = values["mode"].as<int>();
  if(values.count("gene")) {
    auto gmp = construct_gmp(tmp);
    if(!mode) compare_tss(values, gmp, vbed, glist);
    else if(mode==1) merge_gene2bed(values, gmp, vbed, glist);
  } else {
    if(!mode) compare_tss(values, tmp, vbed, glist);
    else if(mode==1) merge_gene2bed(values, tmp, vbed, glist);
  }
  
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

