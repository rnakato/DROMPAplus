/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "dd_command.hpp"
#include "dd_draw.hpp"

namespace {
  void help_global(std::vector<Command> &cmds)
  {
    auto helpmsg = R"(
    ===============

    For the detailed information on the options for each command, use the -h flag along with the command.

    Usage: drompa+ <Command> [options]

    Command:)";
    
    std::cerr << "\n    DROMPA+ v" << VERSION << helpmsg << std::endl;
    for(size_t i=0; i<cmds.size(); ++i) cmds[i].printCommandName();
    std::cerr << std::endl;
    return;
  }
  
  void printVersion()
  {
    std::cerr << "drompa+ version " << VERSION << std::endl;
    exit(0);
  }
  
  int32_t getOpts(std::vector<Command> &cmds, int argc, char* argv[])
  {
    int32_t cmdid(-1);
    
    DEBUGprint("setOpts...");
    
    MyOpt::Opts command("Command");
    MyOpt::Opts genopts("Options");
    
    command.add_options()
      ("command", boost::program_options::value<std::string>(), "command to run");
    genopts.add_options()
      ("version,v", "print version");

    boost::program_options::positional_options_description pd;
    pd.add("command", 1);

    MyOpt::Opts allopts("Options");
    allopts.add(command).add(genopts);

    DEBUGprint("getOpts...");
    
    MyOpt::Variables values;
    
    try {
      // parse first argument only
      boost::program_options::parsed_options parsed = boost::program_options::command_line_parser(2, argv).options(allopts).positional(pd).run();
      store(parsed, values);

      if (values.count("version")) printVersion();
      
      // check command and param
      int32_t on(0);
      std::string cmd = MyOpt::getVal<std::string>(values, "command");
      for (size_t i=0; i<cmds.size(); ++i) {
	if (cmd == cmds[i].getname()) {
	  cmdid = i;
	  cmds[i].SetValue(argc-1, argv+1);
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
      exit(0);
    }
    return cmdid;
  }
  
}

int main(int argc, char* argv[])
{
  auto cmds = generateCommands();

  if (argc == 1) {
    help_global(cmds);
    exit(0);
  }

  int32_t cmdid = getOpts(cmds, argc, argv);
  cmds[cmdid].execute();

  return 0;
}

void MergePdf(DROMPA::Global &p, const std::string &StrAllPdf)
{
  std::string command("cpdf " + StrAllPdf + " -o " + p.getFigFileName());
  std::cout << "Merge PDF files \"" << command << "\"" << std::endl;
  
  int32_t return_code = system(command.c_str());
  if(WEXITSTATUS(return_code)) {
    std::cerr << "Warning: command " << command << "return nonzero status." << std::endl;
  }

  /* rm */
  if (!p.isshowchr()) {
    command = "rm " + StrAllPdf;
    return_code = system(command.c_str());
    if(WEXITSTATUS(return_code)) {
      std::cerr << "Warning: command " << command  << "return nonzero status." << std::endl;
    }
  }
  return;
}

void exec_PCSHARP(DROMPA::Global &p)
{
  std::string StrAllPdf("");
  for(auto &chr: p.gt) {
    if (!p.isincludeYM() && (chr.getname() == "Y" || chr.getname() == "M")) continue;
    if (p.drawregion.getchr() != "" && p.drawregion.getchr() != chr.getname()) continue;
    if (p.drawregion.isRegionLociFile()) {
      int32_t n(0);
      for (auto &m: p.anno.gmp.at(rmchr(chr.getname()))) {
	if(p.drawregion.ExistGeneLociFile(m.second.gname)) ++n;
      }
#ifdef DEBUG
      printf("DEBUG geneloci\n");
      std::cout << "chr" << chr.getname() << ": " << std::flush;
      printf("n=%d\n",0);
#endif
      if(!n) continue;
    }

    std::cout << "chr" << chr.getname() << ": " << std::flush;

    std::vector<bed> regionBed(p.drawregion.getRegionBedChr(chr.getname()));
    if (p.drawregion.isRegionBed() && !regionBed.size()) continue;
 
    Figure fig(p, chr);
    if (fig.Draw(p)) StrAllPdf += p.getFigFileNameChr(chr.getrefname()) + " ";
  }

  MergePdf(p, StrAllPdf);
  return;
}

void exec_GV(DROMPA::Global &p)
{
  p.isGV = true;
  int32_t lenmax(0);
  for (auto &x: p.gt) {
    if (lenmax < x.getlen()) lenmax = x.getlen();
  }
  p.drawparam.width_per_line = lenmax + 1;

  p.drawregion.isRegionOff();
  
  std::string StrAllPdf("");
  for(auto &chr: p.gt) {
    if(!p.isincludeYM() && (chr.getname() == "Y" || chr.getname() == "M")) continue;

    Figure fig(p, chr);
    if (fig.Draw(p)) StrAllPdf += p.getFigFileNameChr(chr.getrefname()) + " ";
  }

  MergePdf(p, StrAllPdf);
  return;
}

template <class T>
void MakeProfile(DROMPA::Global &p)
{
  T profile(p);
  profile.setOutputFilename(p);
  profile.printHead(p);

  for(auto &chr: p.gt) {
    if(!p.isincludeYM() && (chr.getname() == "Y" || chr.getname() == "M")) continue;
    std::cout << "\nchr" << chr.getname() << "..";

    profile.WriteTSV_EachChr(p, chr);

    DEBUGprint("WriteTSV_EachChr done.");
  }

  profile.printNumOfSites();
  profile.MakeFigure(p);
  return;
}

void exec_PROFILE(DROMPA::Global &p)
{  
  if (p.prof.isPtypeTSS() || p.prof.isPtypeTTS()) MakeProfile<ProfileTSS>(p);
  else if (p.prof.isPtypeGene100()) MakeProfile<ProfileGene100>(p);
  else if (p.prof.isPtypeBed())     MakeProfile<ProfileBedSites>(p);

  return;
}
