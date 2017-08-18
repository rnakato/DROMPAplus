/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _UTIL_HPP_
#define _UTIL_HPP_
#include <iostream>
#include <string>
#include <vector>
#include <boost/format.hpp>
#include <boost/program_options.hpp>

std::string rmchr(const std::string &chr);
void isFile(const std::string &);
int32_t isStr(std::string, std::string);


void printList();

template <class Thead, class... Tbody>
void printList(Thead head, Tbody... body)
{
  std::cout << head;
  if(sizeof...(body) > 0) std::cout << '\t';
  printList(body...);
}


template <class T>
void chkminus(const boost::program_options::variables_map &values, std::string x, int32_t lim)
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

template <class T>
void chkrange(const boost::program_options::variables_map &values, std::string x, T s, T e)
{
  if (values.count(x)) {
    T val = values[x].as<T>();
    if(!my_range(val, s, e)) {
      std::cerr << "Error: invalid " << x << ": " << val << std::endl;
      exit(1);
    }
  }
  return;
  }

template <class T>
int32_t getmaxi(std::vector<T> v)
{
  T max(0);
  int32_t maxi(0);
  int32_t size = v.size();
  for(int32_t i=0; i<size; ++i) {
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
  int32_t size = v.size();

  std::vector<double> w(4,0);
  double var=1;
  for(int32_t j=0; j<4; ++j) w[j] = exp(static_cast<double>(-j*j)/2*var*var);
  double r = 1/(w[0] + (w[1]+w[2]+w[3]) *2);

  double m0;
  double m1(v[0]);
  double m2(v[1]);
  double m3(v[2]);
  for(int32_t i=3; i<size-3; ++i) {
    m0 = v[i];
    v[i] = (w[0]*m0 + w[1]*(m1 + v[i+1]) + w[2]*(m2 + v[i+2]) + w[3]*(m3 + v[i+3]))*r;
    m3 = m2;
    m2 = m1;
    m1 = m0;
  }
  return;
}

template <class T>
int32_t findIndex(std::vector<T> array, T value)
{
    auto iter = std::find(array.begin(), array.end(), value);
    size_t index = std::distance(array.begin(), iter);
    if(index == array.size()) index = -1;
    
    return index;
}

#endif /* _UTIL_HPP_ */
