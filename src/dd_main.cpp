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
			   {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::THRE, DrompaCommand::ANNO_PC, DrompaCommand::ANNO_GV, DrompaCommand::DRAW, DrompaCommand::REGION, DrompaCommand::SCALE, DrompaCommand::OVERLAY, DrompaCommand::OTHER}));
    cmds.emplace_back(Command("PC_ENRICH","peak-calling (enrichment ratio)",
			   "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			   exec_PCSHARP,
			   {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::THRE, DrompaCommand::ANNO_PC, DrompaCommand::ANNO_GV, DrompaCommand::DRAW, DrompaCommand::REGION, DrompaCommand::SCALE, DrompaCommand::OTHER}));
    cmds.emplace_back(Command("GV", "global-view visualization",
			   "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			   exec_PCSHARP,
			   {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::ANNO_GV, DrompaCommand::DRAW, DrompaCommand::SCALE, DrompaCommand::OTHER}));
    cmds.emplace_back(Command("PD", "peak density",
			   "-pd <pdfile>,<name> [-pd <pdfile>,<name> ...]",
			   //	   dd_pd,
			   exec_PCSHARP,
			   {DrompaCommand::PD, DrompaCommand::ANNO_GV, DrompaCommand::DRAW, DrompaCommand::SCALE, DrompaCommand::OTHER}));
    cmds.emplace_back(Command("CI", "compare peak-intensity between two samples",
			   "-i <ChIP>,,<name> -i <ChIP>,,<name> -bed <bedfile>",
			   exec_PCSHARP,
			   {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::OTHER}));
    cmds.emplace_back(Command("PROFILE", "make R script of averaged read density",
			   "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			   exec_PCSHARP,
			   {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::PROF, DrompaCommand::OTHER}));
    cmds.emplace_back(Command("HEATMAP", "make heatmap of multiple samples",
			   "-i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]",
			   exec_PCSHARP,
			   {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::PROF, DrompaCommand::OTHER}));
    cmds.emplace_back(Command("CG", "output ChIP-reads in each gene body",
			   "-i <ChIP>,,<name> [-i <ChIP>,,<name> ...]",
			   exec_PCSHARP,
			   {DrompaCommand::CHIP, DrompaCommand::CG, DrompaCommand::OTHER}));
    cmds.emplace_back(Command("TR",      "calculate the travelling ratio (pausing index) for each gene",
			   "-i <ChIP>,,<name> [-i <ChIP>,,<name> ...]",
			   exec_PCSHARP,
			   {DrompaCommand::CHIP, DrompaCommand::PROF, DrompaCommand::OTHER}));
    cmds.emplace_back(Command("GOVERLOOK", "genome-wide overlook of peak positions",
			   "-bed <bedfile>,<name> [-bed <bedfile>,<name> ...]",
			   //dd_overlook,
			   exec_PCSHARP,
			   {DrompaCommand::OTHER}));
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
      for(size_t i=0; i<cmds.size(); ++i) {
	if(cmd == cmds[i].name) {
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


void exec_PCSHARP(DROMPA::Global &p)
{
  printf("drompa\n");

  for(auto chr: p.gt) {
    if(!p.isincludeYM() && (chr.getname() == "Y" || chr.getname() == "M")) continue;
    if(p.drawregion.getchr() != "" && p.drawregion.getchr() != chr.getname()) continue;
    //if(d->drawregion_argv && !d->drawregion->chr[chr].num) return;

    std::cout << "chr" << chr.getname() <<": " << p.drawregion.getchr() << std::flush;

    Figure fig(p, chr);
    // readAnnotation();
    printf("drawing...\n");
    fig.DrawData(p);
    //    sprintf(d->command_mergepdf, "%s%s.pdf ", d->command_mergepdf, prefix);
  }

  return;
}

void Command::InitDump()
{
  std::vector<std::string> str_bool = {"OFF", "ON"};
  std::vector<std::string> str_ratio = {"OFF", "Linear", "Logratio"};
  std::vector<std::string> str_gftype = {"refFlat", "Ensembl", "gtf", "SGD"};
  std::vector<std::string> str_wigfiletype = {"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};
  std::vector<std::string> str_norm  = { "OFF", "TOTALREAD", "NCIS" };
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
	DEBUGprint("INITDUMP:DrompaCommand::CHIP");
	std::cout << boost::format("\nSamples\n");
	for(uint i=0; i<p.samplepair.size(); ++i) {
	  std::cout << (i+1) << ": ";
	  p.samplepair[i].print();
	}
	if (values.count("if")) std::cout << boost::format("Input format: %1%\n") % str_wigfiletype[MyOpt::getVal<int32_t>(values, "if")];
	break;
      }
    case DrompaCommand::NORM:
      {
	DEBUGprint("INITDUMP:DrompaCommand::NORM");
	std::cout << boost::format("   ChIP/Input normalization: %s\n") % str_norm[MyOpt::getVal<int32_t>(values, "norm")];
	if(MyOpt::getVal<int32_t>(values, "sm")) std::cout << boost::format("   smoothing width: %1% bp\n") % MyOpt::getVal<int32_t>(values, "sm");
	break;
      }
    case DrompaCommand::THRE: 
      {
	DEBUGprint("INITDUMP:DrompaCommand::THRE");
	/*	std::cout << boost::format("   Peak intensity threshold: %1$.2f\n")               % values["IPmaxthre"].as<double>();
	std::cout << boost::format("   Enrichment threshold: %1$.2f\n")                   % values["enrichthre"].as<double>();
	std::cout << boost::format("   p-value threshold (internal, -log10): %1$.2e\n")   % values["pthre_internal"].as<double>();
	std::cout << boost::format("   p-value threshold (enrichment, -log10): %1$.2e\n") % values["pthre_enrich"].as<double>();
	std::cout << boost::format("   FDR threshold: %1$.2e\n")                          % values["FDR"].as<double>();*/
  	break;
      }
    case DrompaCommand::ANNO_PC:
      {
	DEBUGprint("INITDUMP:DrompaCommand::ANNO_PC");
	std::cout << boost::format("\nAnnotations:\n");
	if(values.count("gene")) std::cout << boost::format("   Gene file: %1%, Format: %2%\n")
				   % MyOpt::getVal<std::string>(values, "gene")
				   % str_gftype[MyOpt::getVal<int32_t>(values, "gftype")];
	MyOpt::printOpt<std::string>(values, "ars",    "   ARS file");
	MyOpt::printOpt<std::string>(values, "ter",    "   TER file");
	MyOpt::printOpt<std::string>(values, "repeat", "   Repeat file");
	MyOpt::printOpt<std::string>(values, "gc", "   GCcontents file");
	MyOpt::printOpt<std::string>(values, "gd", "   Gene density file");
	/*	if(d->arsfile)     std::cout << boost::format("   ARS file: %1%\n")          % values["ars"].as<std::string>();
	if(d->terfile)     std::cout << boost::format("   TER file: %1%\n")          % values["ter"].as<std::string>();
	if(d->repeat.argv) std::cout << boost::format("   Repeat file: %1%\n")       % values["repeat"].as<std::string>();*/
	MyOpt::printOpt<std::string>(values, "region", "   Region file");
	MyOpt::printVOpt<std::string>(values, "bed", "   Bed file");
	//	if(name != "PROFILE" || name != "HEATMAP") std::cout << boost::format("   name: %1%\n") % d->bed[i]->name;
	break;
      }
    case DrompaCommand::ANNO_GV:
      {
	DEBUGprint("INITDUMP:DrompaCommand::ANNO_GV");
	std::cout << boost::format("\nAnnotations:\n");
	MyOpt::printOpt<std::string>(values, "gc", "   GCcontents file");
	MyOpt::printOpt<std::string>(values, "gd", "   Gene density file");
	/*	
	if(d->GC.argv)     std::cout << boost::format("   GCcontents file: %1%\n")   % values["gc"].as<std::string>();
	if(d->GD.argv)     std::cout << boost::format("   Gene density file: %1%\n") % values["gd"].as<std::string>();*/
	MyOpt::printVOpt<std::string>(values, "inter", "   Interaction file");
	if (values.count("mp")) {
	  std::cout << boost::format("Mappability file directory: %1%\n") % MyOpt::getVal<std::string>(values, "mp");
	  std::cout << boost::format("\tLow mappablitiy threshold: %1%\n") % MyOpt::getVal<double>(values, "mpthre");
	}
	break;
      }
    case DrompaCommand::DRAW:
      {
	DEBUGprint("INITDUMP:DrompaCommand::DRAW");
	std::cout << boost::format("\nFigure parameter:\n");
	std::cout << boost::format("   Display read: ChIP %1%, Input %2%\n") % str_bool[MyOpt::getVal<int32_t>(values, "showctag")] % str_bool[MyOpt::getVal<int32_t>(values, "showitag")];
	std::cout << boost::format("   Display enrichment: %1%\n")           % str_ratio[MyOpt::getVal<int32_t>(values, "showratio")];
	std::cout << boost::format("   Display pvalue (internal): %1%\n")    % str_bool[MyOpt::getVal<int32_t>(values, "showpinter")];
	std::cout << boost::format("   Display pvalue (ChIP/Input): %1%\n")  % str_bool[MyOpt::getVal<int32_t>(values, "showpenrich")];
	std::cout << boost::format("   Y label: %1%\n")                      % str_bool[!values.count("offylabel")];
	std::cout << boost::format("   Y memory: %1%\n")                     % str_bool[!values.count("offymem")];
	break;
      }
    
    case DrompaCommand::REGION:
      {
	DEBUGprint("INITDUMP:DrompaCommand::REGION");
	break;
      }
    case DrompaCommand::SCALE:
      {
	DEBUGprint("INITDUMP:DrompaCommand::SCALE");
	for (auto x: {"scale_tag","scale_ratio","scale_pvalue","scale_tag2","scale_ratio2","scale_pvalue2","bn","ystep"}) chkminus<int>(values, x, 0);
	break;
      }
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
      
    case DrompaCommand::OTHER:
      {
	DEBUGprint("INITDUMP:DrompaCommand::OTHER");
	break;
      }
    }
  }
  
  if(values.count("chr")) std::cout << boost::format("output %1% only.\n") % MyOpt::getVal<int32_t>(values, "chr");

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
	DEBUGprint("Global::setValues::ChIP");
	if (!values.count("input")) PRINTERR("specify --input option.");
	std::vector<std::string> v(MyOpt::getVal<std::vector<std::string>>(values, "input"));
	for(auto x:v) scan_samplestr(x, sample, samplepair, iftype);

	if (values.count("if")) iftype = static_cast<WigType>(MyOpt::getVal<int32_t>(values, "if")); 
	break;
      }
    case DrompaCommand::NORM:
      {
	DEBUGprint("Global::setValues::NORM");
	chkrange<int>(values, "norm", 0, 1);
	break;
      }
    case DrompaCommand::THRE: 
      {
	DEBUGprint("Global::setValues::THRE");
	/*	for (auto x: {"pthre_internal", "pthre_enrich", "qthre", "ipm", "ethre"}) chkminus<int>(values, x, -1);
		chkminus<int>(values, "width4lmd", 0);*/
	break;
      }
    case DrompaCommand::ANNO_PC:
      {
	DEBUGprint("Global::setValues::ANNO_PC");
	chkrange<int>(values, "gftype", 0, 3);
	for (auto x: {"gene", "ars", "ter"}) if (values.count(x)) isFile(MyOpt::getVal<std::string>(values, x));
	break;
      }
    case DrompaCommand::ANNO_GV:
      {
	DEBUGprint("Global::setValues::ANNO_GV");
	/*	chkminus<int>(values, "mpthre", -1);
		for (auto x: {"gcsize", "gdsize"}) chkminus<int>(values, x, 0);*/
	break;
      }
    case DrompaCommand::DRAW:
      {
	DEBUGprint("Global::setValues::DRAW");
	drawparam.setValues(values, samplepair.size());
	break;
      }  
    case DrompaCommand::REGION:
      {
	DEBUGprint("Global::setValues::REGION");
	drawregion.setValues(values);
	break;
      }
    case DrompaCommand::SCALE:
      {
	DEBUGprint("Global::setValues::SCALE");
	scale.setValues(values);
	drawparam.setlineheight(scale.getlineheight());
	break;
      }
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
    case DrompaCommand::OTHER:
      {
	DEBUGprint("Global::setValues::OTHER");
	if (values.count("includeYM")) includeYM = true;
	if (values.count("png"))       ispng = true;
	break;
      }
    }
  }
  return;
}



void DROMPA::Global::setOpts(std::vector<DrompaCommand> &st)
{

  MyOpt::Opts o("Required",100);
  o.add_options()
    ("output,o",  boost::program_options::value<std::string>(), "Output prefix")
    ("gt",        boost::program_options::value<std::string>(), "Genome table")
    ;
  opts.add(o);
  for(auto x: st) {
    switch(x) {
    case DrompaCommand::CHIP:
      {
	boost::program_options::options_description o("Input",100);
	o.add_options()
	  ("input,i",
	   boost::program_options::value<std::vector<std::string>>(),
	   "Specify ChIP data, Input data and name of ChIP sample\n     (separated by ',', values except for 1 can be omitted)\n     1:ChIP   2:Input   3:name   4:peaklist   5:binsize\n     6:scale_tag   7:scale_ratio   8:scale_pvalue\n")
	  ("if",
	   boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::range<int32_t>, _1, 0, static_cast<int>(WigType::WIGTYPENUM) -2, "--if")),
	   "Input file format\n   0: binary (.bin)\n   1: compressed wig (.wig.gz)\n   2: uncompressed wig (.wig)\n   3: bedGraph (.bedGraph)\n   4: bigWig (.bw)")
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::NORM:
      {
	boost::program_options::options_description o("",100);
	o.add_options()
	  ("norm",      boost::program_options::value<int32_t>()->default_value(1),	     "Normalization between ChIP and Input\n      0: not normalize\n      1: with total read number\n      2: with NCIS method\n")
	  ("sm",        boost::program_options::value<int32_t>()->default_value(0),      "Smoothing width") // gausian ??
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::THRE: 
      {
	boost::program_options::options_description o("Threshold",100);
	o.add_options()
	  ("pthre_internal", boost::program_options::value<double>()->default_value(1e-4), "p-value for ChIP internal")
	  ("pthre_enrich",   boost::program_options::value<double>()->default_value(1e-4), "p-value for ChIP/Input enrichment")
	  ("qthre",          boost::program_options::value<double>()->default_value(1),    "FDR")
	  ("ethre,e",        boost::program_options::value<double>()->default_value(2),    "IP/Input fold enrichment")
	  ("ipm",            boost::program_options::value<double>()->default_value(0),    "Read intensity of peak summit")
	  ("nosig", "Omit highlighting peak regions")
	  ("width4lmd", boost::program_options::value<int32_t>()->default_value(100000), "Width for calculating local lambda")
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::ANNO_PC:
      {
	boost::program_options::options_description o("Annotation",100);
	o.add_options()
	  ("gene,g", boost::program_options::value<std::string>(),	  "Gene annotation file")
	  ("gftype", boost::program_options::value<int32_t>()->default_value(1), "Format of gene annotation\n     0: RefFlat (default)\n     1: Ensembl\n     2: GTF (for S. pombe)\n     3: SGD (for S. cerevisiae)\n")
	  ("ars",    boost::program_options::value<std::string>(),	  "ARS list (for yeast)")
	  ("ter",    boost::program_options::value<std::string>(),	  "TER list (for S.cerevisiae)")  
	  ("bed",    boost::program_options::value<std::vector<std::string>>(), "<bedfile>,<label>: Specify bed file and name (<label> can be omited)")
	  ("repeat", boost::program_options::value<std::string>(),	  "Display repeat annotation (RepeatMasker format)") 
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::ANNO_GV:
      {
	boost::program_options::options_description o("Optional data",100);
	o.add_options()
	  ("mp",     boost::program_options::value<std::string>(),  	  "Mappability file")
	  ("mpthre", boost::program_options::value<double>()->default_value(0.3), "Low mappability threshold")
	  ("gap",    boost::program_options::value<std::string>(),	  "Specify gapped regions to be shaded")
	  ("inter",  boost::program_options::value<std::vector<std::string>>(), "<interaction file>,<label>: Specify interaction file and name (<label> can be omited)")  // FDRde iro kaeru
	  ("gc",     boost::program_options::value<std::string>(), 	  "Visualize GC contents graph")
	  ("gcsize", boost::program_options::value<int32_t>()->default_value(100000), "Window size for GC contents")
	  ("gd",     boost::program_options::value<std::string>(), 	  "Visualize gene density (number of genes for each window)")
	  ("gdsize", boost::program_options::value<int32_t>()->default_value(100000), "Window size for gene density")
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::DRAW:   drawparam.setOpts(opts); break;
    case DrompaCommand::REGION: drawregion.setOpts(opts); break;
    case DrompaCommand::SCALE:  scale.setOpts(opts); break;
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
      
    case DrompaCommand::OTHER:
      {
	boost::program_options::options_description o("Others",100);
	o.add_options()
	  ("includeYM", "output peaks of chromosome Y and M")
	  ("rmchr",   "Remove chromosome-separated pdf files")
	  ("png",     "Output with png format (Note: output each page separately)")
	  ("threads,p",
	   boost::program_options::value<int32_t>()->default_value(1)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--thread")),
	   "number of threads to launch")
	  ("help,h", "show help message")
	  ;
	opts.add(o);
	break;
      }
    }
  }

}


/*
void opt::add(std::vector<DrompaCommand> st)
{
  boost::program_options::options_description o("Required",100);
  o.add_options()
    ("output,o",  boost::program_options::value<std::string>(),	 "Output prefix")
    ("gt",        boost::program_options::value<std::string>(),	 "Genome table")
    ;
  opts.add(o);

  for(auto x: st) {
    switch(x) {
    case DrompaCommand::CHIP:
      {
	boost::program_options::options_description o("Input",100);
	o.add_options()
	  ("input,i",   boost::program_options::value<std::vector<std::string>>(), "Specify ChIP data, Input data and name of ChIP sample\n     (separated by ',', values except for 1 can be omitted)\n     1:ChIP   2:Input   3:name   4:peaklist   5:binsize\n     6:scale_tag   7:scale_ratio   8:scale_pvalue\n")
	("if", boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::range<int32_t>, _1, 0, static_cast<int>(WigType::WIGTYPENUM) -2, "--if")),
	 "Input file format\n   0: binary (.bin)\n   1: compressed wig (.wig.gz)\n   2: uncompressed wig (.wig)\n   3: bedGraph (.bedGraph)\n   4: bigWig (.bw)")
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::NORM:
      {
	boost::program_options::options_description o("",100);
	o.add_options()
	  ("norm",      boost::program_options::value<int32_t>()->default_value(1),	     "Normalization between ChIP and Input\n      0: not normalize\n      1: with total read number\n      2: with NCIS method\n")
	  ("sm",        boost::program_options::value<int32_t>()->default_value(0),      "Smoothing width") // gausian ??
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::THRE: 
      {
	boost::program_options::options_description o("Threshold",100);
	o.add_options()
	  ("pthre_internal", boost::program_options::value<double>()->default_value(1e-4), "p-value for ChIP internal")
	  ("pthre_enrich",   boost::program_options::value<double>()->default_value(1e-4), "p-value for ChIP/Input enrichment")
	  ("qthre",          boost::program_options::value<double>()->default_value(1),    "FDR")
	  ("ethre,e",        boost::program_options::value<double>()->default_value(2),    "IP/Input fold enrichment")
	  ("ipm",            boost::program_options::value<double>()->default_value(0),    "Read intensity of peak summit")
	  ("nosig", "Omit highlighting peak regions")
	  ("width4lmd", boost::program_options::value<int32_t>()->default_value(100000), "Width for calculating local lambda")
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::ANNO_PC:
      {
	boost::program_options::options_description o("Annotation",100);
	o.add_options()
	  ("gene,g", boost::program_options::value<std::string>(),	  "Gene annotation file")
	  ("gftype", boost::program_options::value<int32_t>()->default_value(1), "Format of gene annotation\n     0: RefFlat (default)\n     1: Ensembl\n     2: GTF (for S. pombe)\n     3: SGD (for S. cerevisiae)\n")
	  ("ars",    boost::program_options::value<std::string>(),	  "ARS list (for yeast)")
	  ("ter",    boost::program_options::value<std::string>(),	  "TER list (for S.cerevisiae)")  
	  ("bed",    boost::program_options::value<std::vector<std::string>>(), "<bedfile>,<label>: Specify bed file and name (<label> can be omited)")
	  ("repeat", boost::program_options::value<std::string>(),	  "Display repeat annotation (RepeatMasker format)") 
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::ANNO_GV:
      {
	boost::program_options::options_description o("Optional data",100);
	o.add_options()
	  ("mp",     boost::program_options::value<std::string>(),  	  "Mappability file")
	  ("mpthre", boost::program_options::value<double>()->default_value(0.3), "Low mappability threshold")
	  ("gap",    boost::program_options::value<std::string>(),	  "Specify gapped regions to be shaded")
	  ("inter",  boost::program_options::value<std::vector<std::string>>(), "<interaction file>,<label>: Specify interaction file and name (<label> can be omited)")  // FDRde iro kaeru
	  ("gc",     boost::program_options::value<std::string>(), 	  "Visualize GC contents graph")
	  ("gcsize", boost::program_options::value<int32_t>()->default_value(100000), "Window size for GC contents")
	  ("gd",     boost::program_options::value<std::string>(), 	  "Visualize gene density (number of genes for each window)")
	  ("gdsize", boost::program_options::value<int32_t>()->default_value(100000), "Window size for gene density")
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::DRAW:
      {
	opts.add();
	/*	boost::program_options::options_description o("Drawing",100);
	o.add_options()
	  ("showctag",     boost::program_options::value<int32_t>()->default_value(1),    "Display ChIP read lines")
	  ("showitag",     boost::program_options::value<int32_t>()->default_value(0),    "Display Input read lines (0:off 1:all 2:first one)")
	  ("showratio",    boost::program_options::value<int32_t>()->default_value(0),    "Display ChIP/Input ratio (0:off 1:liner scale 2:logscale)")
	  ("showpinter",   boost::program_options::value<int32_t>()->default_value(0),    "Display -log10(p) lines for ChIP internal")
	  ("showpenrich",  boost::program_options::value<int32_t>()->default_value(0),    "Display -log10(p) lines for ChIP/Input enrichment")
	  ("showars",     boost::program_options::value<int32_t>()->default_value(0),     "Display ARS only (do not display genes)")
	  ("ls",          boost::program_options::value<int32_t>()->default_value(1000), "Width for each line (kb)")
	  ("lpp",         boost::program_options::value<int32_t>()->default_value(1),    "Line number per page")
	  ("offymem",     boost::program_options::value<int32_t>()->default_value(0),     "Omit Y memory")
	  ("offylab",     boost::program_options::value<int32_t>()->default_value(0),     "Omit Y label")
	  ("viz",         boost::program_options::value<int32_t>()->default_value(0), "Color of read profile\n     0: normal color\n     1: semitransparent color\n")
	  ;
	  opts.add(o);*/
	
/*	break;
	}
    case DrompaCommand::REGION:
      {
	boost::program_options::options_description o("Region to draw",100);
	o.add_options()
	  ("chr",         boost::program_options::value<int32_t>(),     "Output the specified chromosome only")
	  ("region,r",    boost::program_options::value<std::string>(),  "Specify genomic regions for drawing")
	  ("genefile",    boost::program_options::value<std::string>(),  "Specify gene loci to visualize")  
	  ("len_genefile",boost::program_options::value<int32_t>()->default_value(50000), "extended length for each gene locus")
	  ;
	opts.add(o);
	break;
      }
    
    case DrompaCommand::SCALE:
      {
	boost::program_options::options_description o("Scale for Y axis",100);
	o.add_options()
	  ("scale_tag",    boost::program_options::value<double>(), "Scale for read line")
	  ("scale_ratio",  boost::program_options::value<double>(), "Scale for fold enrichment")
	  ("scale_pvalue", boost::program_options::value<double>(), "Scale for -log10(p)")
	  ("bn",           boost::program_options::value<int32_t>()->default_value(2),     "Number of memories of y-axis")
	  ("ystep",        boost::program_options::value<double>()->default_value(20), "Height of read line")
	  ;
	opts.add(o);
	break;
      }
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
      
    case DrompaCommand::OTHER:
      {
	boost::program_options::options_description o("Others",100);
	o.add_options()
	  ("includeYM", "output peaks of chromosome Y and M")
	  ("rmchr",   "Remove chromosome-separated pdf files")
	  ("png",     "Output with png format (Note: output each page separately)")
	  ("threads,p",
	   boost::program_options::value<int32_t>()->default_value(1)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--thread")),
	   "number of threads to launch")
	  ("help,h", "show help message")
	  ;
	opts.add(o);
	break;
      }
    }
  }
  }*/


/*

void Command::checkParam() {
  for (auto x: {"output", "gt"}) if (!values.count(x)) PRINTERR("specify --" << x << " option.");

  p.getoprefix(MyOpt::getVal<std::string>(values, "output"));
  p.gt = read_genometable(MyOpt::getVal<std::string>(values, "gt"));
	
  for(auto op: vopts) {
    switch(op) {
    case DrompaCommand::CHIP:
      {
	for (auto x: {"input"}) if (!values.count(x)) PRINTERR("specify --" << x << " option.");

	std::vector<std::string> v(MyOpt::getVal<std::vector<std::string>>(values, "input"));
	for(auto x:v) scan_samplestr(x, p);

	if (values.count("if")) p.getWigType(MyOpt::getVal<int32_t>(values, "if")); 
	break;
      }
    case DrompaCommand::NORM:
      {
	chkminus<int>(values, "sm", 0);
	chkrange<int>(values, "norm", 0, 1);
	break;
      }
    case DrompaCommand::THRE: 
      {
	for (auto x: {"pthre_internal", "pthre_enrich", "qthre", "ipm", "ethre"}) chkminus<int>(values, x, -1);
	chkminus<int>(values, "width4lmd", 0);
	break;
      }
    case DrompaCommand::ANNO_PC:
      {
	chkrange<int>(values, "gftype", 0, 3);
	for (auto x: {"gene", "ars", "ter"}) if (values.count(x)) isFile(MyOpt::getVal<std::string>(values, x));
	break;
      }
    case DrompaCommand::ANNO_GV:
      {
	chkminus<int>(values, "mpthre", -1);
	for (auto x: {"gcsize", "gdsize"}) chkminus<int>(values, x, 0);
	break;
      }
    case DrompaCommand::DRAW:
      {
	for (auto x: {"ls", "lpp"}) chkminus<int>(values, x, 0);
	break;
      }
    case DrompaCommand::REGION:
      {
	for (auto x: {"region", "genefile"}) if (values.count(x)) isFile(MyOpt::getVal<std::string>(values, x));
	chkminus<int>(values, "len_genefile", -1);
	break;
      }
    case DrompaCommand::SCALE:
      {
	for (auto x: {"scale_tag","scale_ratio","scale_pvalue","bn","ystep"}) chkminus<int>(values, x, 0);
	break;
      }
    case DrompaCommand::OVERLAY:
      {
	for (auto x: {"scale_tag2","scale_ratio2","scale_pvalue2"}) chkminus<int>(values, x, 0);
	break;
      }
    case DrompaCommand::CG: 
      {
	for (auto x: {"cgthre"}) chkminus<int>(values, x, -1);
	break;
      }
    case DrompaCommand::TR: 
      {
	for (auto x: {"tssthre"}) chkminus<int>(values, x, -1);
	break;
      }
    case DrompaCommand::PD:
      {
	for (auto x: {"pd"}) if (!values.count(x)) PRINTERR("specify --" << x << " option.");
	for (auto x: {"pdsize"}) chkminus<int>(values, x, 0);
	
	std::vector<std::string> v(MyOpt::getVal<std::vector<std::string>>(values, "pd"));
	for(auto &x: v) p.pd.emplace_back(scan_pdstr(x));
	break;
      }
    case DrompaCommand::PROF:
      {
	chkrange<int>(values, "ptype", 0, 4);
	chkrange<int>(values, "stype", 0, 1);
	chkrange<int>(values, "ntype", 0, 1);
	for (auto x: {"cw", "maxval", "hmsort"}) chkminus<int>(values, x, 0);
	break;
      }
      
    case DrompaCommand::OTHER:
      {
	if (values.count("includeYM")) p.includeYM = true;
	if (values.count("png"))       p.setpng();
	break;
      }
    }
    return;
  }
}

*/
