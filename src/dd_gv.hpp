/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_GV_H_
#define _DD_GV_H_

#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <boost/format.hpp>
#include "SSP/common/inline.hpp"
#include "SSP/common/seq.hpp"
#include "SSP/common/util.hpp"
#include "dd_class.hpp"
#include "version.hpp"

enum class DrompaCommand {CHIP, NORM, THRE, ANNO_PC, ANNO_GV, DRAW, REGION, SCALE, CG, PD, TR, PROF, OVERLAY, OTHER};

class opt {
public:
  boost::program_options::options_description opts;
  opt(const std::string str): opts(str) {}
  void add(std::vector<DrompaCommand> st);
};

  class Command {
    opt opts;
    std::string desc;
    std::string requiredstr;
    std::vector<DrompaCommand> vopts;
    boost::program_options::variables_map values;
    //    std::function<void(boost::program_options::variables_map &, Param &)> func;
  
    std::vector<chrsize> gt;
    std::unordered_map<std::string, SampleFile> sample;
    std::vector<SamplePair> samplepair;
    std::vector<pdSample> pd;

  public:
    std::string name;

    Command(std::string n, std::string d, std::string r,
	    //	    std::function<void(boost::program_options::variables_map &, Param &)> _func,
	    std::vector<DrompaCommand> v): opts("Options"), desc(d), requiredstr(r), vopts(v),
					   //func(_func),
					   name(n) {
      opts.add(v);
    };
    void print() const {
      std::cout << std::setw(8) << " " << std::left << std::setw(12) << name
		<< std::left << std::setw(40) << desc << std::endl;
    }
    void printhelp() const {
      std::cout << boost::format("%1%:  %2%\n") % name % desc;
      std::cout << boost::format("Usage: drompa %1% [options] -o <output> --gt <genometable> %2%\n\n") % name % requiredstr;
      std::cout << opts.opts << std::endl;
    }
    void checkParam();
    void InitDump();
    void SetValue(int argc, char* argv[]) {
      if (argc ==1) {
	printhelp();
	exit(0);
      }
      try {
	store(parse_command_line(argc, argv, opts.opts), values);
      
	if (values.count("help")) {
	  printhelp();
	  exit(0);
	}

	notify(values);
	checkParam();
	InitDump();

      } catch (std::exception &e) {
	std::cout << e.what() << std::endl;
      }
    }
    void execute() {
      //      func(values, p);
    }
  };


#endif /* _DD_GV_H_ */
