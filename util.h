/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _UTIL_H_
#define _UTIL_H_
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/format.hpp>

template <class T>
void chkminus(const boost::program_options::variables_map &values, std::string x, int lim)
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
  void printr(std::ofstream &out, T a, S b)
{
  double r = b ? a*100/(double)b: 0;
  out << boost::format("%1% (%2$.1f%%)\t") % a % r;
};

template <class T>
int getmaxi(std::vector<T> v)
{
  T max(0);
  int maxi(0);
  int size = v.size();
  for(int i=0; i<size; ++i) {
    if(max < v[i]) {
      max = v[i];
      maxi = i;
    }
  }
  return maxi;
};

#endif /* _UTIL_H_ */
