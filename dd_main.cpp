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


//             PC_BROAD    peak-calling (for broad mode)
//             FRIP        accumulate read counts in bed regions specified
//             3DMAP       accumulate read counts in bed regions specified

vector<Command> cmds({
    {"PC_SHARP", "peak-calling (for sharp mode)",
	{OPTREQ, OPTIO, OPTTHRE, OPTANNO, OPTDRAW, OPTSCALE, OPTOTHER}},
    {"PC_ENRICH","peak-calling (enrichment ratio)",
	{OPTREQ, OPTIO, OPTTHRE, OPTANNO, OPTDRAW, OPTSCALE, OPTOTHER}},
      {"GV",
	    "global-view visualization",
	      {OPTREQ, OPTIO, OPTTHRE, OPTANNO, OPTDRAW, OPTSCALE, OPTOTHER}},
	  {"PD",
	      "peak density",
		{OPTREQ, OPTIO, OPTTHRE, OPTANNO, OPTDRAW, OPTSCALE, OPTOTHER}},
	    {"CI",
		"compare peak-intensity between two samples",
		  {OPTREQ, OPTIO, OPTOTHER}},
	      {"CG",
		  "output ChIP-reads in each gene body",
		    {OPTREQ, OPTIO, OPTCG, OPTOTHER}},
		{"GOVERLOOK",
		    "genome-wide overlook of peak positions",
		      {OPTREQ, OPTIO, OPTOTHER}},
	    {"PROFILE",   "make R script of averaged read density",
		{OPTREQ, OPTIO, OPTPROF, OPTOTHER}},
	      {"HEATMAP",   "make heatmap of multiple samples", {OPTREQ, OPTIO, OPTPROF, OPTOTHER}},
		{"TR",        "calculate the travelling ratio (pausing index) for each gene", {OPTREQ, OPTIO, OPTPROF, OPTOTHER}},
});


void help_global() {
    auto helpmsg = R"(
    ===============

    For the detailed information on the options for each command, use the -h flag along with the command.

    Usage: drompa_draw <Command> [options]

    Command:)";

    cerr << "\n    DROMPA v" << VERSION << helpmsg << endl;
    for(size_t i=0; i<cmds.size(); ++i) cmds[i].print();
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

    // check command and param
    int on(0);
    string cmd = values["command"].as<string>();
    for(size_t i=0; i<cmds.size(); ++i) {
      if(cmd == cmds[i].name) {
	cmds[i].getOpts(argc-1, argv+1);
	on++;
      }
    }
    
    if (!on) {
      cerr << "  Invalid command: " << values["command"].as<string>() << endl;
      help_global();
      exit(0);
    }

  } catch (exception &e) {
    cout << e.what() << endl;
  }
  
  return 0;
}
