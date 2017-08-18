/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <time.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "ssp_gv.hpp"
#include "ParseMapfile.hpp"

#define VERSION "1.0.1"

namespace {
  void printVersion()
  {
    std::cerr << "SSP version " << VERSION << std::endl;
    exit(0);
  }
  
  void help_global()
  {
    auto helpmsg = R"(
===============

Usage: ssp [option] -i <inputfile> -o <output> --gt <genome_table>)";
    
    std::cerr << "\nSSP v" << VERSION << helpmsg << std::endl;
    return;
  }
  
  void estimateFragLength(SSP::Global &p)
  {
    DEBUGprint("estimateFragLength...");
    if(p.genome.dflen.isnomodel()) return; // p.genome.isPaired()
    
    std::string head(p.getprefix());
    
    clock_t t1,t2;
    t1 = clock();
    strShiftProfile(p.sspst, p.genome, head, "jaccard"); 
    t2 = clock();
    std::cout << "Jaccard Bit: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
    makeFCSProfile(p.fcsst, p.genome, head, "fcs");
    p.outputSSPstats();
    
    clock_t t3 = clock();
    std::cout << "Fragment variability: " << static_cast<double>(t3 - t2) / CLOCKS_PER_SEC << "sec.\n";
    
    if(p.sspst.DoExjac()) {
      t1 = clock();
      strShiftProfile(p.sspst, p.genome, head, "exjaccard");
      t2 = clock();
      std::cout << "Jaccard Vec: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
    }

    if(p.sspst.DoHd()) {
      t1 = clock();
      strShiftProfile(p.sspst, p.genome, head, "hdp");
      t2 = clock();
      std::cout << "Hamming: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
    }
  
    if(p.sspst.DoCc()) {
      t1 = clock();
      strShiftProfile(p.sspst, p.genome, head, "ccp");
      t2 = clock();    
      std::cout << "ccp: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
    }
    
    DEBUGprint("estimateFragLength done.");
    return;
  }

  void init_dump(const MyOpt::Variables &values, SSP::Global &ssp){
    std::vector<std::string> str_wigfiletype = {"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};
 
    std::cout << boost::format("\n======================================\n");
    std::cout << boost::format("SSP version %1%\n\n") % VERSION;

    MyOpt::dumpIO(values);
    MyOpt::dumpGenomeTable(values);
  
    MyOpt::dumpFragmentLengthDist(values);
    MyOpt::dumpPair(values);
    MyOpt::dumpLibComp(values);
    ssp.sspst.dump();
    
    MyOpt::dumpOther(values);
    printf("======================================\n");
    return;
  }
  
  void getOpts(SSP::Global &ssp, int argc, char* argv[])
  {
    DEBUGprint("setOpts...");

    MyOpt::Opts allopts("Options");
    MyOpt::setOptIO(allopts, "sspout");
    MyOpt::setOptPair(allopts);
    ssp.setOpts(allopts);
    MyOpt::setOptOther(allopts);

    DEBUGprint("getOpts...");
    
    MyOpt::Variables values;
    
    try {
      boost::program_options::parsed_options parsed = parse_command_line(argc, argv, allopts);
      store(parsed, values);
    }
    catch(const boost::program_options::error_with_option_name& e) {
      std::cout << e.what() << std::endl;
      exit(0);
    }
    if (values.count("version")) printVersion();
    if (argc ==1) {
      help_global();
      std::cerr << "Use --help option for more information on the other options\n\n";
      exit(0);
    }
    if (values.count("help")) {
      help_global();
      std::cout << "\n" << allopts << std::endl;
      exit(0);
    }
    std::vector<std::string> opts = {"input", "output", "gt"};
    for (auto x: opts) {
      if (!values.count(x)) PRINTERR("specify --" << x << " option.");
    }
      
    try {
      notify(values);
      ssp.setValues(values);
    
      boost::filesystem::path dir(MyOpt::getVal<std::string>(values, "odir"));
      boost::filesystem::create_directory(dir);
    
      init_dump(values, ssp);
    } catch(const boost::bad_any_cast& e) {
      std::cout << e.what() << std::endl;
      exit(0);
    }
    
    DEBUGprint("getOpts done.");
    return;
  }
  
}

int main(int argc, char* argv[])
{
  SSP::Global p;
  getOpts(p, argc, argv);

  read_mapfile(p.genome);
  p.complexity.checkRedundantReads(p.genome);

  estimateFragLength(p);

#ifdef DEBUG
  p.genome.printReadstats();
#endif

  return 0;
}
