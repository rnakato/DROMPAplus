/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <fstream>
#include <unordered_map>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <ext/stdio_filebuf.h>
#include "dd_gv.hpp"
#include "dd_opt.hpp"
#include "WigStats.hpp"
#include "SSP/common/ReadAnnotation.hpp"
#include "SSP/src/SeqStats.hpp"
#include "SSP/common/inline.hpp"

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
			 {DrompaOpt::CHIP, DrompaOpt::NORM, DrompaOpt::THRE, DrompaOpt::ANNO_PC, DrompaOpt::ANNO_GV, DrompaOpt::DRAW, DrompaOpt::REGION, DrompaOpt::SCALE, DrompaOpt::OVERLAY, DrompaOpt::OTHER}));
  cmds.push_back(Command("PC_ENRICH","peak-calling (enrichment ratio)",
			 "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			 drompa,
			 {DrompaOpt::CHIP, DrompaOpt::NORM, DrompaOpt::THRE, DrompaOpt::ANNO_PC, DrompaOpt::ANNO_GV, DrompaOpt::DRAW, DrompaOpt::REGION, DrompaOpt::SCALE, DrompaOpt::OTHER}));
  cmds.push_back(Command("GV", "global-view visualization",
			 "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			 drompa,
			 {DrompaOpt::CHIP, DrompaOpt::NORM, DrompaOpt::ANNO_GV, DrompaOpt::DRAW, DrompaOpt::SCALE, DrompaOpt::OTHER}));
  cmds.push_back(Command("PD", "peak density",
			 "-pd <pdfile>,<name> [-pd <pdfile>,<name> ...]",
			 dd_pd,
			 {DrompaOpt::PD, DrompaOpt::ANNO_GV, DrompaOpt::DRAW, DrompaOpt::SCALE, DrompaOpt::OTHER}));
  cmds.push_back(Command("CI", "compare peak-intensity between two samples",
			 "-i <ChIP>,,<name> -i <ChIP>,,<name> -bed <bedfile>",
			 drompa,
			 {DrompaOpt::CHIP, DrompaOpt::NORM, DrompaOpt::OTHER}));
  cmds.push_back(Command("PROFILE", "make R script of averaged read density",
			 "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			 drompa,
			 {DrompaOpt::CHIP, DrompaOpt::NORM, DrompaOpt::PROF, DrompaOpt::OTHER}));
  cmds.push_back(Command("HEATMAP", "make heatmap of multiple samples",
			 "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			 drompa,
			 {DrompaOpt::CHIP, DrompaOpt::NORM, DrompaOpt::PROF, DrompaOpt::OTHER}));
  cmds.push_back(Command("CG", "output ChIP-reads in each gene body",
			 "-i <ChIP>,,<name> [-i <ChIP>,,<name> ...]",
			 drompa,
			 {DrompaOpt::CHIP, DrompaOpt::CG, DrompaOpt::OTHER}));
  cmds.push_back(Command("TR",      "calculate the travelling ratio (pausing index) for each gene",
			 "-i <ChIP>,,<name> [-i <ChIP>,,<name> ...]",
			 drompa,
			 {DrompaOpt::CHIP, DrompaOpt::PROF, DrompaOpt::OTHER}));
  cmds.push_back(Command("GOVERLOOK", "genome-wide overlook of peak positions",
			 "-bed <bedfile>,<name> [-bed <bedfile>,<name> ...]",
			 dd_overlook,
			 {DrompaOpt::OTHER}));
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
    std::string cmd = MyOpt::getVal<std::string>(values, "command");
    for(size_t i=0; i<cmds.size(); ++i) {
      if(cmd == cmds[i].name) {
	cmds[i].execute(argc-1, argv+1);
	on++;
      }
    }

    if (!on) {
      std::cerr << "  Invalid command: " << MyOpt::getVal<std::string>(values, "command") << std::endl;
      help_global(cmds);
      exit(0);
    }

  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}

template <class T>
void readWig(T &in, WigArray &array, std::string filename, std::string chrname, int binsize)
{
  std::string head("chrom="+ chrname +"\tspan=");
  int on(0);

  std::string lineStr;
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || !lineStr.find("track")) continue;
    if(on && isStr(lineStr, "chrom=")) break;
    if(isStr(lineStr, head)) {
      on=1;
      continue;
    }
    if(!on) continue;
    std::vector<std::string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    int i = (stoi(v[0])-1)/binsize;
    array.setval(i, stol(v[1]));
  }

  // array.printArray();
  return;
}

void readBinary(WigArray &array, const std::string &filename, const int32_t nbin)
{
  static int nbinsum(0);
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (!in) PRINTERR("cannot open " << filename);

  in.seekg(nbinsum * sizeof(int32_t));

  array.readBinary(in, nbin);
  // array.printArray();

  nbinsum += nbin;
  return;
}

WigArray read_wigdata(variables_map &values, std::unordered_map<std::string, SampleFile>::iterator itr, chrsize &chr)
{
  std::cout << chr.getname() << std::endl;
  int binsize(itr->second.getbinsize());
  int nbin(chr.getlen()/binsize +1);
  std::string filename = itr->first;
  
  WigArray array(nbin, 0);
  WigType iftype(itr->second.getiftype());

  if (iftype == WigType::UNCOMPRESSWIG) {
    std::ifstream in(filename);
    if (!in) PRINTERR("cannot open " << filename);
    readWig(in, array, filename, chr.getname(), binsize);
  } else if (iftype == WigType::COMPRESSWIG) {
    std::string command = "zcat " + filename;
    FILE *fp = popen(command.c_str(), "r");
    __gnu_cxx::stdio_filebuf<char> *p_fb = new __gnu_cxx::stdio_filebuf<char>(fp, std::ios_base::in);
    std::istream in(static_cast<std::streambuf *>(p_fb));
    readWig(in, array, filename, chr.getname(), binsize);
  } else if (iftype == WigType::BEDGRAPH) {
    //    outputBedGraph(values, p, filename);
  } else if (iftype == WigType::BINARY) {
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
      auto wigarray = read_wigdata(values, itr, chr);
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
