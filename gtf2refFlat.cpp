#include <string>
#include <iostream>
#include <algorithm>
#include "cmdline.h"
#include "readdata.h"

using namespace std;

cmdline::parser argv_init(int argc, char* argv[])
{
  cmdline::parser p;
  p.add<string>("gtf", 'g', "gtf file", true, "");
  p.add("name", 'n', "output name instead of id");
  p.add("help", 'h', "print this message");

  if (argc==1 || !p.parse(argc, argv) || p.exist("help")) {
    if (argc==1 || p.exist("help")) cout << "Parser for gtf file." << endl;
    cout << p.error_full() << p.usage();
    exit(1);
  }
  return p;
}

int main(int argc, char* argv[])
{
  cmdline::parser p = argv_init(argc, argv);

  auto tmp = parseGtf(p.get<string>("gtf"), p.exist("name"));  // hash for transcripts
  auto gmp = construct_gmp(tmp);                 // hash for genes

  //printMap(tmp);
  //printMap(gmp);
  printRefFlat(tmp);

  return 0;
}
