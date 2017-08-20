/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_GV_H_
#define _DD_GV_H_

#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <boost/format.hpp>
#include "dd_class.hpp"
#include "version.hpp"
#include "SSP/common/inline.hpp"
#include "SSP/common/seq.hpp"
#include "SSP/common/util.hpp"
#include "SSP/common/BoostOptions.hpp"

class Command {
  std::string desc;
  std::string requiredstr;
  std::vector<DrompaCommand> vopts;
  boost::program_options::variables_map values;
  std::function<void(DROMPA::Global &p)> func;

  DROMPA::Global p;

public:
  std::string name;

  Command(std::string n, std::string d, std::string r,
	  std::function<void(DROMPA::Global &p)> _func,
	  std::vector<DrompaCommand> v): desc(d), requiredstr(r), vopts(v),
					 func(_func), name(n)
  {
    p.setOpts(v);
  };
  void printCommandName() const {
    std::cout << std::setw(8) << " " << std::left << std::setw(12) << name
	      << std::left << std::setw(40) << desc << std::endl;
  }
  void printhelp() const {
    std::cout << boost::format("%1%:  %2%\n") % name % desc;
    std::cout << boost::format("Usage: drompa %1% [options] -o <output> --gt <genometable> %2%\n\n") % name % requiredstr;
    std::cout << p.opts << std::endl;
  }
  void InitDump();
  void SetValue(int argc, char* argv[]) {
    if (argc ==1) {
      printhelp();
      exit(0);
    }
    try {
      store(parse_command_line(argc, argv, p.opts), values);
      
      if (values.count("help")) {
	printhelp();
	exit(0);
      }

      notify(values);
      p.setValues(vopts, values);
      InitDump();

    } catch (std::exception &e) {
      std::cout << e.what() << std::endl;
    }
  }

  void execute() {
    func(p);
  }
};


#endif /* _DD_GV_H_ */
