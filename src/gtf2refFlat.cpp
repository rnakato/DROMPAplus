#include <string>
#include <iostream>
#include <algorithm>
#include <boost/program_options.hpp>
#include "SSP/src/readdata.h"

using namespace std;
using namespace boost::program_options;

variables_map argv_init(int argc, char* argv[])
{
  options_description allopts("Options");
  allopts.add_options()
    ("gtf,g", value<string>(), "Gene file")
    ("name,n", "Output name instead of id")
    ("unique,u", "Only output one transcript per one gene (default: all transcripts)")
    ("help,h", "Print this message")
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
    for (auto x: {"gtf"}) {
      if (!values.count(x)) {
	cerr << "specify --" << x << " option." << endl;
	exit(0);
      }
    }
    notify(values);

  } catch (exception &e) {
    cout << e.what() << endl;
    exit(0);
  }
  return values;
}

int main(int argc, char* argv[])
{
  variables_map values = argv_init(argc, argv);

  auto tmp = parseGtf(values["gtf"].as<string>(), values.count("name")); // hash for transcripts
  auto gmp = construct_gmp(tmp);                 // hash for genes

  //printMap(tmp);
  //printMap(gmp);
    
  if (values.count("unique")) printRefFlat(gmp); 
  else printRefFlat(tmp);

  return 0;
}
