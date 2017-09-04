/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <boost/format.hpp>
#include "dd_gv.hpp"
#include "dd_readfile.hpp"
#include "dd_draw.hpp"
#include "WigStats.hpp"
#include "SSP/common/BoostOptions.hpp"
#include "SSP/common/ReadAnnotation.hpp"
#include "SSP/src/SeqStats.hpp"
#include "SSP/common/inline.hpp"

//             PC_BROAD    peak-calling (for broad mode)
//             FRIP        accumulate read counts in bed regions specified
//             3DMAP       accumulate read counts in bed regions specified

void exec_PCSHARP(DROMPA::Global &p);
void exec_GV(DROMPA::Global &p);


void help_global(std::vector<Command> &cmds) {
    auto helpmsg = R"(
    ===============

    For the detailed information on the options for each command, use the -h flag along with the command.

    Usage: drompa_draw <Command> [options]

    Command:)";

    std::cerr << "\n    DROMPA v" << VERSION << helpmsg << std::endl;
    for(size_t i=0; i<cmds.size(); ++i) cmds[i].printCommandName();
      std::cerr << std::endl;
    return;
}

void printVersion()
{
  std::cerr << "drompa version " << VERSION << std::endl;
  exit(0);
}

namespace {
  
  std::vector<Command> generateCommands()
  {
    std::vector<Command> cmds;
    cmds.emplace_back(Command("PC_SHARP", "peak-calling (for sharp mode)",
			   "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			   exec_PCSHARP,
			      {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::THRE, DrompaCommand::ANNO_PC, DrompaCommand::ANNO_GV, DrompaCommand::DRAW, DrompaCommand::REGION, DrompaCommand::OVERLAY, DrompaCommand::OTHER},
			      CommandParamSet(5, 1, 0, 30, 4, 5)));
    cmds.emplace_back(Command("PC_ENRICH","peak-calling (enrichment ratio)",
			   "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			   exec_PCSHARP,
			   {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::THRE, DrompaCommand::ANNO_PC, DrompaCommand::ANNO_GV, DrompaCommand::DRAW, DrompaCommand::REGION, DrompaCommand::OTHER},
			      CommandParamSet(5, 0, 1, 30, 4, 5)));
    cmds.emplace_back(Command("GV", "global-view visualization",
			   "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			   exec_GV,
			   {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::ANNO_GV, DrompaCommand::DRAW, DrompaCommand::OTHER},
			      CommandParamSet(0, 0, 1, 2000, 2, 10)));
    cmds.emplace_back(Command("PD", "peak density",
			   "-pd <pdfile>,<name> [-pd <pdfile>,<name> ...]",
			   //	   dd_pd,
			   exec_PCSHARP,
			   {DrompaCommand::PD, DrompaCommand::ANNO_GV, DrompaCommand::DRAW, DrompaCommand::OTHER},
			      CommandParamSet(0, 0, 0, 0, 0, 0)));
    cmds.emplace_back(Command("CI", "compare peak-intensity between two samples",
			   "-i <ChIP>,,<name> -i <ChIP>,,<name> -bed <bedfile>",
			   exec_PCSHARP,
			   {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::OTHER},
			      CommandParamSet(0, 0, 0, 0, 0, 0)));
    cmds.emplace_back(Command("PROFILE", "make R script of averaged read density",
			   "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			   exec_PCSHARP,
			   {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::PROF, DrompaCommand::OTHER},
			      CommandParamSet(0, 0, 0, 0, 0, 0)));
    cmds.emplace_back(Command("HEATMAP", "make heatmap of multiple samples",
			   "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			   exec_PCSHARP,
			   {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::PROF, DrompaCommand::OTHER},
			      CommandParamSet(0, 0, 0, 30, 4, 5)));
    cmds.emplace_back(Command("CG", "output ChIP-reads in each gene body",
			   "-i <ChIP>,,<name> [-i <ChIP>,,<name> ...]",
			   exec_PCSHARP,
			   {DrompaCommand::CHIP, DrompaCommand::CG, DrompaCommand::OTHER},
			      CommandParamSet(0, 0, 0, 0, 0, 0)));
    cmds.emplace_back(Command("TR",      "calculate the travelling ratio (pausing index) for each gene",
			   "-i <ChIP>,,<name> [-i <ChIP>,,<name> ...]",
			   exec_PCSHARP,
			   {DrompaCommand::CHIP, DrompaCommand::PROF, DrompaCommand::OTHER},
			      CommandParamSet(0, 0, 0, 0, 0, 0)));
    cmds.emplace_back(Command("GOVERLOOK", "genome-wide overlook of peak positions",
			   "-bed <bedfile>,<name> [-bed <bedfile>,<name> ...]",
			   //dd_overlook,
			   exec_PCSHARP,
			   {DrompaCommand::OTHER},
			      CommandParamSet(5, 0, 0, 0, 0, 0)));
    return cmds;
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
  printf("Drawing...\n");

  std::string StrAllPdf("");
  for(auto &chr: p.gt) {
    if(!p.isincludeYM() && (chr.getname() == "Y" || chr.getname() == "M")) continue;
    if(p.drawregion.getchr() != "" && p.drawregion.getchr() != chr.getname()) continue;

    Figure fig(p, chr);
    if (fig.Draw(p, chr)) StrAllPdf += p.getFigFileNameChr(chr.getrefname()) + " ";
  }

  MergePdf(p, StrAllPdf);
  return;
}

void exec_GV(DROMPA::Global &p)
{
  printf("Drawing...\n");
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
    if (fig.Draw(p, chr)) StrAllPdf += p.getFigFileNameChr(chr.getrefname()) + " ";
  }

  MergePdf(p, StrAllPdf);
  return;
}

void Command::InitDump()
{
  std::vector<std::string> str_bool = {"OFF", "ON"};
  std::vector<std::string> str_stype = { "ChIP read", "Enrichment ratio", "Enrichment P-value" };
  std::vector<std::string> str_ptype = { "NONE", "TSS", "TTS", "GENE100", "SPECIFIEDSITES" };
  std::vector<std::string> str_ntype = { "WHOLE GENOME", "TARGET REGIONS ONLY" };

  std::cout << boost::format("\n======================================\n");
  std::cout << boost::format("drompa version %1%: %2%\n\n") % VERSION % name;
  std::cout << boost::format("output prefix: %1%\n")     % MyOpt::getVal<std::string>(values, "output");
  std::cout << boost::format("Genome-table file: %1%\n") % MyOpt::getVal<std::string>(values, "gt");

  for(auto x: vopts) {
    switch(x) {
    case DrompaCommand::CHIP:
      {
	std::vector<std::string> str_wigfiletype = {"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};
	DEBUGprint("INITDUMP:DrompaCommand::CHIP");
	std::cout << boost::format("\nSamples\n");
	for(uint i=0; i<p.samplepair.size(); ++i) {
	  std::cout << (i+1) << ": ";
	  p.samplepair[i].print();
	}
	if (MyOpt::getVal<int32_t>(values, "if") < static_cast<int32_t>(WigType::NONE)) std::cout << boost::format("Input format: %1%\n") % str_wigfiletype[MyOpt::getVal<int32_t>(values, "if")];
	break;
      }
    case DrompaCommand::NORM: p.InitDumpNorm(); break;
    case DrompaCommand::THRE: p.thre.InitDump(); break;
    case DrompaCommand::ANNO_PC: p.anno.InitDumpPC(values); break;
    case DrompaCommand::ANNO_GV: p.anno.InitDumpGV(values); break;
    case DrompaCommand::DRAW: p.drawparam.InitDump(); break;
    case DrompaCommand::REGION: p.drawregion.InitDump(values); break;
    case DrompaCommand::OVERLAY:
      {
	DEBUGprint("INITDUMP:DrompaCommand::OVERLAY");
	if(values.count("ioverlay")) { 
	  std::cout << boost::format("\nOverlayed samples\n");
	  for(uint i=0; i<p.samplepair.size(); ++i) {
	    std::cout << (i+1) << ": ";
	    p.samplepair[i].print();
	  }
	}
	break;
      }
    case DrompaCommand::CG: 
      {
	DEBUGprint("INITDUMP:DrompaCommand::CG");
	for (auto x: {"cgthre"}) chkminus<int32_t>(values, x, -1);
	break;
      }
    case DrompaCommand::TR: 
      {
	DEBUGprint("INITDUMP:DrompaCommand::TR");
	for (auto x: {"tssthre"}) chkminus<int32_t>(values, x, -1);
	break;
      }
    case DrompaCommand::PD:
      {
	DEBUGprint("INITDUMP:DrompaCommand::PD");
	std::cout << boost::format("\nSamples\n");
	for(uint i=0; i<p.pd.size(); ++i) {
	  std::cout << boost::format("   IP%1%: %2%\tname: %3%\n") % (i+1) % p.pd[i].argv % p.pd[i].name;
	}
	break;
      }
    case DrompaCommand::PROF:
      {
	DEBUGprint("INITDUMP:DrompaCommand::PROF");
	std::cout << boost::format("   show type: %1$\n")             % str_stype[MyOpt::getVal<int32_t>(values, "stype")];
	std::cout << boost::format("   profile type: %1$\n")          % str_ptype[MyOpt::getVal<int32_t>(values, "ptype")];
	std::cout << boost::format("   profile normalization: %1$\n") % str_ntype[MyOpt::getVal<int32_t>(values, "ntype")];
	break;
      }
      
    case DrompaCommand::OTHER: p.InitDumpOther(); break;
    }
  }
  
  printf("======================================\n");
  return;
}


void DROMPA::Global::setValues(const std::vector<DrompaCommand> &vopts, const MyOpt::Variables &values)
{
  for (auto x: {"output", "gt"}) if (!values.count(x)) PRINTERR("specify --" << x << " option.");

  oprefix = MyOpt::getVal<std::string>(values, "output");
  gt = read_genometable(MyOpt::getVal<std::string>(values, "gt"));
	
  for(auto op: vopts) {
    switch(op) {
    case DrompaCommand::CHIP:
      {
	DEBUGprint("ChIP setValues...");
	try {
	  if (!values.count("input")) PRINTERR("specify --input option.");
	  if (values.count("if")) iftype = static_cast<WigType>(MyOpt::getVal<int32_t>(values, "if"));
	  std::vector<std::string> v(MyOpt::getVal<std::vector<std::string>>(values, "input"));
	  for(auto x:v) scan_samplestr(x, gt, sample, samplepair, iftype);

	} catch(const boost::bad_any_cast& e) {
	  std::cout << e.what() << std::endl;
	  exit(0);
	}
	DEBUGprint("ChIP setValues done.");
  
	break;
      }
    case DrompaCommand::NORM: setValuesNorm(values); break;
    case DrompaCommand::THRE: thre.setValues(values); break;
    case DrompaCommand::ANNO_PC: anno.setValuesPC(values); break;
    case DrompaCommand::ANNO_GV: anno.setValuesGV(values); break;
    case DrompaCommand::DRAW: drawparam.setValues(values, samplepair.size()); break;
    case DrompaCommand::REGION: drawregion.setValues(values); break;
    case DrompaCommand::OVERLAY:
      {
	DEBUGprint("Global::setValues::OVERLAY");
	for (auto x: {"scale_tag2","scale_ratio2","scale_pvalue2"}) chkminus<int>(values, x, 0);
	break;
      }
    case DrompaCommand::CG: 
      {
	DEBUGprint("Global::setValues::CG");
	for (auto x: {"cgthre"}) chkminus<int>(values, x, -1);
	break;
      }
    case DrompaCommand::TR: 
      {
	DEBUGprint("Global::setValues::TR");
	for (auto x: {"tssthre"}) chkminus<int>(values, x, -1);
	break;
      }
    case DrompaCommand::PD:
      {
	DEBUGprint("Global::setValues::PD");
	for (auto x: {"pd"}) if (!values.count(x)) PRINTERR("specify --" << x << " option.");
	for (auto x: {"pdsize"}) chkminus<int>(values, x, 0);
	
	std::vector<std::string> v(MyOpt::getVal<std::vector<std::string>>(values, "pd"));
	for(auto &x: v) pd.emplace_back(scan_pdstr(x));
	break;
      }
    case DrompaCommand::PROF:
      {
	DEBUGprint("Global::setValues::PROF");
	chkrange<int>(values, "ptype", 0, 4);
	chkrange<int>(values, "stype", 0, 1);
	chkrange<int>(values, "ntype", 0, 1);
	for (auto x: {"cw", "maxval", "hmsort"}) chkminus<int>(values, x, 0);
	break;
      }
    case DrompaCommand::OTHER: setValuesOther(values); break;
    }
  }
  return;
}


void DROMPA::Global::setOpts(const std::vector<DrompaCommand> &st, const CommandParamSet &cps)
{

  MyOpt::Opts o("Required",100);
  o.add_options()
    ("output,o",  boost::program_options::value<std::string>(), "Output prefix")
    ("gt",        boost::program_options::value<std::string>(), "Genome table")
    ;
  opts.add(o);
  for(auto &x: st) {
    switch(x) {
    case DrompaCommand::CHIP:
      {
	boost::program_options::options_description o("Input",100);
	o.add_options()
	  ("input,i",
	   boost::program_options::value<std::vector<std::string>>(),
	   "Specify ChIP data, Input data and name of ChIP sample\n     (separated by ',', values except for 1 can be omitted)\n     1:ChIP   2:Input   3:name   4:peaklist   5:binsize\n     6:scale_tag   7:scale_ratio   8:scale_pvalue\n")
	  ("if",
	   boost::program_options::value<int32_t>()->default_value(static_cast<int32_t>(WigType::NONE))->notifier(boost::bind(&MyOpt::range<int32_t>, _1, 0, static_cast<int32_t>(WigType::WIGTYPENUM) -1, "--if")),
	   "Input file format\n   0: binary (.bin)\n   1: compressed wig (.wig.gz)\n   2: uncompressed wig (.wig)\n   3: bedGraph (.bedGraph)\n   4: bigWig (.bw)")
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::NORM:    setOptsNorm(opts, cps.sm); break;
    case DrompaCommand::THRE:    thre.setOpts(opts); break;
    case DrompaCommand::ANNO_PC: anno.setOptsPC(opts); break;
    case DrompaCommand::ANNO_GV: anno.setOptsGV(opts); break;
    case DrompaCommand::DRAW:    drawparam.setOpts(opts, cps); break;
    case DrompaCommand::REGION:  drawregion.setOpts(opts); break;
    case DrompaCommand::OVERLAY:
      {
	boost::program_options::options_description o("For overlay",100);
	o.add_options()
	  ("ioverlay",  boost::program_options::value<std::vector<std::string>>(), "Input file")
	  ("scale_tag2",   boost::program_options::value<double>(), "Scale for read line")
	  ("scale_ratio2", boost::program_options::value<double>(), "Scale for fold enrichment")
	  ("scale_pvalue2",boost::program_options::value<double>(), "Scale for -log10(p)")
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::CG: 
      {
	boost::program_options::options_description o("CG",100);
	o.add_options()
	  ("cgthre",    boost::program_options::value<double>(), "Minimum threshold per kbp")
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::TR: 
      {
	boost::program_options::options_description o("TR",100);
	o.add_options()
	  ("tssthre",    boost::program_options::value<double>(), "")
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::PD:
      {
	boost::program_options::options_description o("PD",100);
	o.add_options()
	  ("pd",   boost::program_options::value<std::vector<std::string>>(), "Peak density file and name\n(separated by ',' <name> can be omited)")
	  ("prop",   boost::program_options::value<double>(),  "scale_tag")
	  ("pdsize", boost::program_options::value<int32_t>()->default_value(100000), "windowsize for peak density")
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::PROF:
      {
	boost::program_options::options_description o("PROFILE AND HEATMAP",100);
	o.add_options()
	  ("ptype",   boost::program_options::value<int32_t>(),  "Region type: 1; around TSS, 2; around TES, 3; divide gene into 100 subregions 4; around peak sites")
	  ("stype",   boost::program_options::value<int32_t>(),  "Show type: 0; ChIP read (default) 1; ChIP/Input enrichment")
	  ("ntype",   boost::program_options::value<int32_t>(),  "Normalization type: 0; total read 1; target regions only")
	  ("cw",      boost::program_options::value<double>()->default_value(2500), "width from the center")
	  ("maxval",   boost::program_options::value<double>(),  "Upper limit for heatmap")
	  ("offse",  "Omit the standard error in profile")
	  ("hmsort",   boost::program_options::value<int32_t>()->default_value(1),  "Column number for sorting sites")
	  ("sortgbody",  "Sort sites by read number of gene body (default: TSS)")
	  ("pdetail",  "")
	  ;
	opts.add(o);
	break;
      }
      
    case DrompaCommand::OTHER: setOptsOther(opts); break;
    }
  }

}
