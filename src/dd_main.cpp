/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <ext/stdio_filebuf.h>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include "SSP/src/readdata.h"
#include "SSP/src/macro.h"
#include "dd_gv.h"
#include "dd_opt.h"

using namespace boost::program_options;

//             PC_BROAD    peak-calling (for broad mode)
//             FRIP        accumulate read counts in bed regions specified
//             3DMAP       accumulate read counts in bed regions specified

void drompa(variables_map &, Param &);
void dd_pd(variables_map &, Param &);
void dd_overlook(variables_map &, Param &);

void help_global(std::vector<Command> &cmds) {
    auto helpmsg = R"(
    ===============

    For the detailed information on the options for each command, use the -h flag along with the command.

    Usage: drompa_draw <Command> [options]

    Command:)";

    std::cerr << "\n    DROMPA v" << VERSION << helpmsg << std::endl;
    for(size_t i=0; i<cmds.size(); ++i) cmds[i].print();
      std::cerr << std::endl;
    return;
}

void printVersion()
{
  std::cerr << "drompa version " << VERSION << std::endl;
  exit(0);
}

std::vector<Command> generateCommands()
{
  std::vector<Command> cmds;
  cmds.push_back(Command("PC_SHARP", "peak-calling (for sharp mode)",
			 "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			 drompa,
			 {OPTCHIP, OPTNORM, OPTTHRE, OPTANNO_PC, OPTANNO_GV, OPTDRAW, OPTREGION, OPTSCALE, OPTOVERLAY, OPTOTHER}));
  cmds.push_back(Command("PC_ENRICH","peak-calling (enrichment ratio)",
			 "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			 drompa,
			 {OPTCHIP, OPTNORM, OPTTHRE, OPTANNO_PC, OPTANNO_GV, OPTDRAW, OPTREGION, OPTSCALE, OPTOTHER}));
  cmds.push_back(Command("GV", "global-view visualization",
			 "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			 drompa,
			 {OPTCHIP, OPTNORM, OPTANNO_GV, OPTDRAW, OPTSCALE, OPTOTHER}));
  cmds.push_back(Command("PD", "peak density",
			 "-pd <pdfile>,<name> [-pd <pdfile>,<name> ...]",
			 dd_pd,
			 {OPTPD, OPTANNO_GV, OPTDRAW, OPTSCALE, OPTOTHER}));
  cmds.push_back(Command("CI", "compare peak-intensity between two samples",
			 "-i <ChIP>,,<name> -i <ChIP>,,<name> -bed <bedfile>",
			 drompa,
			 {OPTCHIP, OPTNORM, OPTOTHER}));
  cmds.push_back(Command("PROFILE", "make R script of averaged read density",
			 "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			 drompa,
			 {OPTCHIP, OPTNORM, OPTPROF, OPTOTHER}));
  cmds.push_back(Command("HEATMAP", "make heatmap of multiple samples",
			 "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			 drompa,
			 {OPTCHIP, OPTNORM, OPTPROF, OPTOTHER}));
  cmds.push_back(Command("CG", "output ChIP-reads in each gene body",
			 "-i <ChIP>,,<name> [-i <ChIP>,,<name> ...]",
			 drompa,
			 {OPTCHIP, OPTCG, OPTOTHER}));
  cmds.push_back(Command("TR",      "calculate the travelling ratio (pausing index) for each gene",
			 "-i <ChIP>,,<name> [-i <ChIP>,,<name> ...]",
			 drompa,
			 {OPTCHIP, OPTPROF, OPTOTHER}));
  cmds.push_back(Command("GOVERLOOK", "genome-wide overlook of peak positions",
			 "-bed <bedfile>,<name> [-bed <bedfile>,<name> ...]",
			 dd_overlook,
			 {OPTOTHER}));
  return cmds;
}

int main(int argc, char* argv[])
{
  auto cmds = generateCommands();

  if (argc ==1) {
    help_global(cmds);
    exit(0);
  }

  options_description command("Command");
  options_description genopts("Options");
  
  command.add_options()
    ("command", value<std::string>(), "command to run");
  genopts.add_options()
    ("version,v", "print version");

  positional_options_description pd;
  pd.add("command", 1);

  options_description allopts("Options");
  allopts.add(command).add(genopts);

  variables_map values;

  try {
    // parse first argument only
    parsed_options parsed = command_line_parser(2, argv).options(allopts).positional(pd).run();
    store(parsed, values);

    if (values.count("version")) printVersion();

    // check command and param
    int on(0);
    std::string cmd = values["command"].as<std::string>();
    for(size_t i=0; i<cmds.size(); ++i) {
      if(cmd == cmds[i].name) {
	cmds[i].execute(argc-1, argv+1);
	on++;
      }
    }

    if (!on) {
      std::cerr << "  Invalid command: " << values["command"].as<std::string>() << std::endl;
      help_global(cmds);
      exit(0);
    }

  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}

template <class T>
void readWig(T &in, std::vector<int> &array, std::string filename, std::string chrname, int binsize)
{
  std::string head("chrom="+ chrname +"\tspan=");
  int on(0);

  std::string lineStr;
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || !lineStr.find("track")) continue;
    if(on && lineStr.find("chrom=")!= std::string::npos) break;
    if(lineStr.find(head)!= std::string::npos) {
      on=1;
      continue;
    }
    if(!on) continue;
    std::vector<std::string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    int i = (stoi(v[0])-1)/binsize;
    array[i] = VALUE2WIGARRAY(stol(v[1]));
    if(array[i]) std::cout << i*50 << "\t" << WIGARRAY2VALUE(array[i]) << std::endl;
  }
  return;
}

void readBinary(std::vector<int> &array, std::string filename, int nbin)
{
  static int nbinsum(0);
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (!in) PRINTERR("cannot open " << filename);

  in.seekg(nbinsum*sizeof(int));  
  for(int i=0; i<nbin; ++i) {
    in.read((char *)&array[i], sizeof(int));
    //    if(array[i]) std::cout << i*50 << "\t" << WIGARRAY2VALUE(array[i]) << std::endl;
  }

  nbinsum += nbin;
  return;
}

std::vector<int> read_wigdata(variables_map &values, std::unordered_map<std::string, SampleFile>::iterator itr, chrsize &chr)
{
  std::cout << chr.getname() << std::endl;
  int binsize(itr->second.getbinsize());
  int nbin(chr.getlen()/binsize +1);
  std::vector<int> array(nbin, 0);
  int iftype(itr->second.getiftype());
  std::string filename = itr->first;

  if (iftype==TYPE_UNCOMPRESSWIG) {
    std::ifstream in(filename);
    if (!in) PRINTERR("cannot open " << filename);
    readWig(in, array, filename, chr.getname(), binsize);
  } else if (iftype==TYPE_COMPRESSWIG) {
    std::string command = "zcat " + filename;
    FILE *fp = popen(command.c_str(), "r");
    __gnu_cxx::stdio_filebuf<char> *p_fb = new __gnu_cxx::stdio_filebuf<char>(fp, std::ios_base::in);
    std::istream in(static_cast<std::streambuf *>(p_fb));
    readWig(in, array, filename, chr.getname(), binsize);
  } else if (iftype==TYPE_BEDGRAPH) {
    //    outputBedGraph(values, p, filename);
  } else if (iftype==TYPE_BINARY) {
    readBinary(array, filename, nbin);
  }

  //  if(p->smoothing) smooth_tags(&(s->data), p->smoothing, values["binsize"].as<int>(), chr.nbin);

  return array;
}

void drompa(variables_map &values, Param &p)
{
  printf("drompa\n");

  for(auto chr:p.gt) {
    for(auto itr = p.sample.begin(); itr != p.sample.end(); ++itr) {
      read_wigdata(values, itr, chr);
    }
  }

  return;
}

void dd_pd(variables_map &values, Param &p)
{
  printf("dd_pd\n");
  return;
}

void dd_overlook(variables_map &values, Param &p)
{
  printf("dd_overlook\n");
  return;
}