/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include "BoostOptions.hpp"

namespace MyOpt {
  void setOptIO(Opts &allopts, const std::string &odir_default)
  {
    Opts opt("Input/Output",100);
    opt.add_options()
      ("input,i", boost::program_options::value<std::string>(),
       "Mapping file. Multiple files are allowed (separated by ',')")
      ("output,o",boost::program_options::value<std::string>(),
       "Prefix of output files")
      ("odir",    boost::program_options::value<std::string>()->default_value(odir_default),
       "output directory name")
      ("ftype,f", boost::program_options::value<std::string>(),
       "{SAM|BAM|BOWTIE|TAGALIGN}: format of input file\nTAGALIGN could be gzip'ed (extension: tagAlign.gz)")
      ;
    allopts.add(opt);
  }
  void setOptPair(Opts &allopts)
  {
    Opts opt("For paired-end",100);
    opt.add_options()
      ("pair", "add when the input file is paired-end")
      ("maxins",
       boost::program_options::value<int32_t>()->default_value(500)->notifier(boost::bind(&over<int32_t>, _1, 1, "--maxins")),
       "maximum fragment length")
      ;
    allopts.add(opt);
  }
  void setOptOther(Opts &allopts)
  {
    Opts opt("Others",100);
    opt.add_options()
      ("threads,p",
       boost::program_options::value<int32_t>()->default_value(1)->notifier(boost::bind(&over<int32_t>, _1, 1, "--thread")),
       "number of threads to launch")
      ("version,v", "print version")
      ("help,h", "show help message")
      ;
    allopts.add(opt);
  }

  void dumpIO(const Variables &values)
  {
    std::cout << boost::format("Input file: %1%\n") % values["input"].as<std::string>();
    if(values.count("ftype")) std::cout << boost::format("\tFormat: %1%\n") % values["ftype"].as<std::string>();
    std::cout << boost::format("Output file: %1%/%2%\n") % values["odir"].as<std::string>() % values["output"].as<std::string>();
  }
  
  void dumpGenomeTable(const Variables &values)
  {
    std::cout << boost::format("Genome-table file: %1%\n") % values["gt"].as<std::string>();
    if(values.count("mptable"))
      std::cout << boost::format("genome-table file for mappable bases: %1%\n") % values["mptable"].as<std::string>();
  }
 
  void dumpLibComp(const Variables &values)
  {
    if (!values.count("nofilter")) {
      std::cout << "PCR bias filtering: ON" << std::endl;
      if (values["thre_pb"].as<int32_t>()) std::cout << "PCR bias threshold: > " << values["thre_pb"].as<int32_t>() << std::endl;
      std::cout << boost::format("\t%1% reads used for library complexity\n") % values["ncmp"].as<int64_t>();
    } else {
      std::cout << "PCR bias filtering: OFF" << std::endl;
    }
  }

  void dumpFragmentLengthDist(const Variables &values)
  {
    if (!values.count("pair")) {
      std::cout << "Single-end mode: ";
      std::cout << "fragment length will be estimated by strand-shift profile" << std::endl;
      if (values.count("nomodel")) std::cout << boost::format("Predefined fragment length: %1%\n") % values["flen"].as<int32_t>();
    }
  }


  void dumpPair(const Variables &values)
  {
    if (values.count("pair")) {
      std::cout << "Paired-end mode: ";
      std::cout << boost::format("Maximum fragment length: %1%\n") % values["maxins"].as<int32_t>();
    }
  }

  void dumpOther(const Variables &values)
  {
    std::cout << boost::format("\nNumber of threads: %1%\n") % values["threads"].as<int32_t>();
  }
}
