/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */

#include "dd_command.hpp"

// CommandParamSet
//  int32_t sm;
//  int32_t showctag;
//  int32_t showratio;
//  double scaletag;
//  double scaleratio;
//  double scalepvalue;
//  double thre_pinter;
//  double thre_penrich;
//  double thre_ethre;
//  double thre_ipm;

std::vector<Command> generateCommands()
{
  std::vector<Command> cmds;
  cmds.emplace_back(Command("PC_SHARP", "Visualization: for sharp mark parameter set",
			    "-i <ChIP>,<Input>,<label> [-i <ChIP>,<Input>,<label> ...]",
			    exec_PCSHARP,
			    {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::THRE, DrompaCommand::ANNO_PC, DrompaCommand::ANNO_GV, DrompaCommand::DRAW, DrompaCommand::REGION, DrompaCommand::OTHER},
			    CommandParamSet(5, 1, 0, 30, 3, 5,
					    5, 4, 0, 0)
			    ));
  cmds.emplace_back(Command("PC_BROAD", "Visualization: for broad mark parameter set",
			    "-i <ChIP>,<Input>,<label> [-i <ChIP>,<Input>,<label> ...]",
			    exec_PCSHARP,
			    {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::THRE, DrompaCommand::ANNO_PC, DrompaCommand::ANNO_GV, DrompaCommand::DRAW, DrompaCommand::REGION, DrompaCommand::OTHER},
			    CommandParamSet(20, 1, 0, 30, 3, 5,
					    4, 3, 0, 0)
			    ));
  cmds.emplace_back(Command("PC_ENRICH","Visualization: ChIP/Input enrichment",
			    "-i <ChIP>,<Input>,<label> [-i <ChIP>,<Input>,<label> ...]",
			    exec_PCSHARP,
			    {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::THRE, DrompaCommand::ANNO_PC, DrompaCommand::ANNO_GV, DrompaCommand::DRAW, DrompaCommand::REGION, DrompaCommand::OTHER},
			    CommandParamSet(5, 0, 1, 30, 3, 5,
					    0, 0, 2, 5)
			    ));
  cmds.emplace_back(Command("GV", "Visualization: global-view enrichment",
			    "-i <ChIP>,<Input>,<label> [-i <ChIP>,<Input>,<label> ...]",
			    exec_GV,
			    {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::ANNO_GV, DrompaCommand::DRAW, DrompaCommand::OTHER},
			    CommandParamSet(0, 0, 1, 2000, 1, 10,
					    0, 0, 2, 0)
			    ));
  cmds.emplace_back(Command("PROFILE", "Visualize averaged read density in specified regions",
			    "-i <ChIP>,<Input>,<label> [-i <ChIP>,<Input>,<label> ...]",
			    exec_PROFILE,
			    {DrompaCommand::CHIP, DrompaCommand::ANNO_PC, DrompaCommand::NORM, DrompaCommand::PROF, DrompaCommand::OTHER},
			    CommandParamSet(5, 0, 0, 0, 0, 0,
					    0, 0, 0, 0)
			    ));
  cmds.emplace_back(Command("MULTICI", "Generate matrix of averaged read density in specified BED sites",
			    "-i <ChIP>,<Input>,<label> [-i <ChIP>,<Input>,<label> ...]",
			    exec_MULTICI,
			    {DrompaCommand::CHIP, DrompaCommand::ANNO_PC, DrompaCommand::NORM, DrompaCommand::PROF, DrompaCommand::OTHER},
			    CommandParamSet(5, 0, 0, 0, 0, 0,
					    0, 0, 0, 0)
			    ));
  cmds.emplace_back(Command("GENWIG", "Generate wig data of enrichment or p-value distribution",
			    "-i <ChIP>,<Input>,<label> [-i <ChIP>,<Input>,<label> ...]",
			    exec_GENWIG,
			    {DrompaCommand::CHIP, DrompaCommand::GENWIG, DrompaCommand::NORM},
			    CommandParamSet(5, 0, 0, 0, 0, 0,
					    0, 0, 0, 0)
			    ));
/*  cmds.emplace_back(Command("CI", "compare peak-intensity between two samples",
			    "-i <ChIP>,,<label> -i <ChIP>,,<label> -bed <bedfile>",
			    exec_PCSHARP,
			    {DrompaCommand::CHIP, DrompaCommand::NORM, DrompaCommand::OTHER},
			    CommandParamSet(0, 0, 0, 0, 0, 0, false)));*/
/*  cmds.emplace_back(Command("CG", "output ChIP-reads in each gene body",
			    "-i <ChIP>,,<label> [-i <ChIP>,,<label> ...]",
			    exec_PCSHARP,
			    {DrompaCommand::CHIP, DrompaCommand::CG, DrompaCommand::OTHER},
			    CommandParamSet(0, 0, 0, 0, 0, 0, false)));*/
/*  cmds.emplace_back(Command("TR",      "calculate the travelling ratio (pausing index) for each gene",
			    "-i <ChIP>,,<label> [-i <ChIP>,,<label> ...]",
			    exec_PCSHARP,
			    {DrompaCommand::CHIP, DrompaCommand::PROF, DrompaCommand::OTHER},
			    CommandParamSet(0, 0, 0, 0, 0, 0, false)));*/
  return cmds;
}

void Command::InitDump()
{
  std::vector<std::string> str_bool = {"OFF", "ON"};

  std::cout << boost::format("\n======================================\n");
  std::cout << boost::format("drompa+ version %1%: %2%\n\n") % VERSION % name;
  std::cout << boost::format("output prefix: %1%\n")     % MyOpt::getVal<std::string>(values, "output");
  std::cout << boost::format("Genome-table file: %1%\n") % MyOpt::getVal<std::string>(values, "gt");

  for (auto x: vopts) {
    switch(x) {
    case DrompaCommand::CHIP: p.InitDumpChIP(); break;
    case DrompaCommand::NORM: p.InitDumpNorm(); break;
    case DrompaCommand::THRE: p.thre.InitDump(); break;
    case DrompaCommand::ANNO_PC:
      if (p.drawparam.isshowpdf()) p.anno.InitDumpPC(values);
      break;
    case DrompaCommand::ANNO_GV:
      if (p.drawparam.isshowpdf()) p.anno.InitDumpGV(values);
      break;
    case DrompaCommand::DRAW:
      if (p.drawparam.isshowpdf()) p.drawparam.InitDump();
      break;
    case DrompaCommand::REGION: p.drawregion.InitDump(values); break;
    case DrompaCommand::GENWIG:
      {
	DEBUGprint("INITDUMP:DrompaCommand::GENWIG");
	std::cout << boost::format("Output format %1%\n")
	  % p.genwig_getOutputFileTypeStr();
	break;
      }
    case DrompaCommand::CG:
      {
	DEBUGprint("INITDUMP:DrompaCommand::CG");
	//	for (auto x: {"cgthre"}) chkminus<int32_t>(values, x, -1);
	break;
      }
    case DrompaCommand::TR:
      {
	DEBUGprint("INITDUMP:DrompaCommand::TR");
	//	for (auto x: {"tssthre"}) chkminus<int32_t>(values, x, -1);
	break;
      }
    case DrompaCommand::PROF: p.prof.InitDump(); break;
    case DrompaCommand::OTHER: p.InitDumpOther(); break;
    }
  }

  printf("======================================\n");
  return;
}
