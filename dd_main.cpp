#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include "common.h"
#include "warn.h"
#include "macro.h"
#include "dd_gv.h"
#include "dd_opt.h"

using namespace std;

vector<vector<string>> cmds({
    {"PC_SHARP", "peak-calling (for sharp mode)"},
    {"PC_ENRICH", "peak-calling (enrichment ratio)"},
    {"GV", "global-view visualization"},
    {"PD", "peak density"},
    {"CI", "compare peak-intensity between two samples"},
    {"CG", "output ChIP-reads in each gene body"},
    {"GOVERLOOK", "genome-wide overlook of peak positions"},
    {"PROFILE", "make R script of averaged read density"},
    {"HEATMAP", "make heatmap of multiple samples"},
    {"TR", "calculate the travelling ratio (pausing index) for each gene"},
});

unordered_map<string, function<int(int, char*[])>> cmdmp({
    {"PC_SHARP", func_PCSHARP},
    {"PC_ENRICH", func_ENRICH},
    {"GV", func_GV},
    {"PD", func_PD},
    {"CI", func_CI},
    {"CG", func_CG},
    {"GOVERLOOK", func_GOVERLOOK},
    {"PROFILE", func_PROFILE},
    {"HEATMAP", func_HEATMAP},
    {"TR", func_TR},
});  

void help_global() {
    auto helpmsg = R"(
    ===============

    For the detailed information on the options for each command, use the -h flag along with the command.

    Usage: drompa_draw <Command> [options]

    Command:)";

    cerr << "\n    DROMPA v" << VERSION << helpmsg << endl;
    for(size_t i=0; i<cmds.size(); ++i) {
      cout << setw(8) << " " << left << setw(12) << cmds[i][0]
	   << left << setw(40) << cmds[i][1] << endl;
    }
    cerr << endl;
    return;
}

void printVersion()
{
  cerr << "drompa version " << VERSION << endl;
  exit(0);
}

int main(int argc, char* argv[])
{
  if (argc ==1) {
    help_global();
    exit(0);
  }

  options_description command("Command");
  options_description genopts("Options");
  
  command.add_options()
    ("command", value<string>(), "command to run");
  genopts.add_options()
    ("version,v", "print version")
    ;

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

    auto cmd = cmdmp.find(values["command"].as<string>());
    if (cmd == cmdmp.end()) {
      cerr << "   Invalid command: " << values["command"].as<string>() << endl;
      help_global();
      exit(0);
    } else {
      cmd->second(argc-1, argv+1);
    }

  } catch (exception &e) {
    cout << e.what() << endl;
  }
  
  return 0;
}

variables_map getOpts(int argc, char* argv[], const opt &opts)
{
  if (argc ==1) {
    BPRINT("%1%:  %2%\n") % cmds[2][0] % cmds[2][1];
    BPRINT("Usage: drompa %1% [options] -p <output> -gt <genometable> -i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]\n\n") % argv[0]; 
    cout << opts.opts << endl;
    exit(0);
  }
  
  variables_map values;
  try {
    store(parse_command_line(argc, argv, opts.opts), values);
    
    if (values.count("help")) {
      BPRINT("%1%:  %2%\n") % cmds[2][0] % cmds[2][1];
      BPRINT("Usage: drompa %1% [options] -p <output> -gt <genometable> -i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]\n\n") % argv[0]; 
      cout << opts.opts << endl;
      exit(0);
    }
    
    notify(values);
    cout << values["gt"].as<string>() << endl;
  } catch (exception &e) {
    cout << e.what() << endl;
  }
  
  return values;
}

int func_PCSHARP(int argc, char* argv[])
{
  
  opt allopts("Options");
  allopts.add({OPTREQ, OPTIO, OPTTHRE, OPTANNO, OPTDRAW, OPTSCALE, OPTOTHER});

  variables_map values = getOpts(argc, argv, allopts);

  return 0;
}

int func_ENRICH(int argc, char* argv[])
{
  opt allopts("Options");
  allopts.add({OPTREQ, OPTIO, OPTTHRE, OPTANNO, OPTDRAW, OPTSCALE, OPTOTHER});

  variables_map values = getOpts(argc, argv, allopts);
  return 0;
}
int func_GV(int argc, char* argv[])
{
  opt allopts("Options");
  allopts.add({OPTREQ, OPTIO, OPTTHRE, OPTANNO, OPTDRAW, OPTSCALE, OPTOTHER});

  variables_map values = getOpts(argc, argv, allopts);
  return 0;
}
int func_PD(int argc, char* argv[])
{
  opt allopts("Options");
  allopts.add({OPTREQ, OPTIO, OPTPD, OPTANNO, OPTDRAW, OPTSCALE, OPTOTHER});

  variables_map values = getOpts(argc, argv, allopts);

  return 0;
}
int func_CI(int argc, char* argv[])
{
  opt allopts("Options");
  allopts.add({OPTREQ, OPTIO, OPTOTHER});

  variables_map values = getOpts(argc, argv, allopts);
  return 0;
}
int func_CG(int argc, char* argv[])
{
  opt allopts("Options");
  allopts.add({OPTREQ, OPTIO, OPTCG, OPTOTHER});

  variables_map values = getOpts(argc, argv, allopts);

  return 0;
}
int func_GOVERLOOK(int argc, char* argv[])
{
  opt allopts("Options");
  allopts.add({OPTREQ, OPTIO, OPTOTHER});

  variables_map values = getOpts(argc, argv, allopts);
  return 0;
}
int func_PROFILE(int argc, char* argv[])
{
  opt allopts("Options");
  allopts.add({OPTREQ, OPTIO, OPTPROF, OPTOTHER});

  variables_map values = getOpts(argc, argv, allopts);
  return 0;
}

int func_HEATMAP(int argc, char* argv[])
{
  opt allopts("Options");
  allopts.add({OPTREQ, OPTIO, OPTPROF, OPTOTHER});

  variables_map values = getOpts(argc, argv, allopts);
  return 0;
}

int func_TR(int argc, char* argv[])
{
  opt allopts("Options");
  allopts.add({OPTREQ, OPTIO, OPTTR, OPTOTHER});

  variables_map values = getOpts(argc, argv, allopts);
  return 0;
}
