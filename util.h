/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _UTIL_H_
#define _UTIL_H_
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/format.hpp>


std::string rmchr(const std::string &chr);
void addmp(std::map<int, double> &mpto, const std::map<int, double> &mpfrom, double w=1);

template <class T>
void chkminus(const boost::program_options::variables_map &values, std::string x, int lim)
{
  if (values.count(x)) {
    T val = values[x].as<T>();
    if(val <= lim) {
      std::cerr << "Error: invalid " << x << ": " << val << std::endl;
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

template <class T>
void GaussianSmoothing(std::vector<T> &v)
{
  int size = v.size();

  std::vector<double> w(4,0);
  double var=1;
  for(int j=0; j<4; ++j) w[j] = exp((double)(-j*j)/2*var*var);
  double r = 1/(w[0] + (w[1]+w[2]+w[3]) *2);

  double m0;
  double m1(v[0]);
  double m2(v[1]);
  double m3(v[2]);
  for(int i=3; i<size-3; ++i) {
    m0 = v[i];
    v[i] = (w[0]*m0 + w[1]*(m1 + v[i+1]) + w[2]*(m2 + v[i+2]) + w[3]*(m3 + v[i+3]))*r;
    m3 = m2;
    m2 = m1;
    m1 = m0;
  }
  return;
}

#endif /* _UTIL_H_ */
