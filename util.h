/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _UTIL_H_
#define _UTIL_H_
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

using namespace std;
using namespace boost::program_options;

#define PRINTERR(...)   do{ cerr << "Error: " << __VA_ARGS__ << endl; exit(1); }while(0)

template <class T>
void chkminus(const variables_map &values, string x, int lim)
{
  if (values.count(x)) {
    T val = values[x].as<T>();
    if(val <= lim) {
      cerr << "Error: invalid " << x << ": " << val << endl;
      exit(1);
    }
  }
  return;
}

void isFile(string);

#endif /* _UTIL_H_ */
