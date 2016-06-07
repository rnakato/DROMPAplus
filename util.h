/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _UTIL_H_
#define _UTIL_H_
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/format.hpp>

using namespace std;
using namespace boost::program_options;

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

template <class T, class S>
void printr(ofstream &out, T a, S b)
{
  double r = b ? a*100/(double)b: 0;
  out << boost::format("%1% (%2$.1f%%)\t") % a % r;
};

#endif /* _UTIL_H_ */
