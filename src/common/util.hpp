/* Copyright(c) Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _UTIL_H_
#define _UTIL_H_

#include <string>
#include <vector>
#include <complex>
#include <iostream>

std::string basename(const std::string &path);

template <class T>
void AverageSmoothing(std::vector<T> &v, const int32_t nsmooth)
{
    if (static_cast<int32_t>(v.size()) <= nsmooth) {
        std::cerr << "Warning: the array size (" << std::to_string(v.size()) 
                  << ") is smaller than smoothing width (--sm " << std::to_string(nsmooth) 
                  << "). Smoothing is omitted." << std::endl;
        return;
    }

  std::vector<double> temparray;
  copy(v.begin(), v.end(), std::back_inserter(temparray));

  for (size_t i=0; i<v.size(); ++i) {
    int32_t count(0);
    double ave(0);
    for (int32_t j=-nsmooth; j<=nsmooth; ++j) {
      if (i+j>= 0 && i+j < v.size()) {
        ave += temparray[i+j];
        ++count;
      }
    }
    if (count) v[i] = ave/count; // else v[i] = v[i];
  }

  return;
}

template <class T>
void GaussianSmoothing(std::vector<T> &v, const int32_t nsmooth)
{
    if (static_cast<int32_t>(v.size()) <= nsmooth) {
        std::cerr << "Warning: the array size (" << std::to_string(v.size()) 
                  << ") is smaller than smoothing width (--sm " << std::to_string(nsmooth) 
                  << "). Smoothing is omitted." << std::endl;
        return;
    }

  std::vector<double> w(nsmooth + 1, 0);
  double var(1);

  double sum(0);
  for (int32_t i=0; i<=nsmooth; ++i) {
    w[i] = exp(static_cast<double>(-i*i)/2*var*var);
    sum += w[i];
  }
  double r(1/(sum*2 - w[0]));

  std::vector<double> m(nsmooth+1,0);
  for (int32_t i=0; i<=nsmooth; ++i) m[i] = v[nsmooth-i];

  //  std::cout << "GS Weight: ";
  //for (int32_t i=0; i<=nsmooth; ++i)  std::cout << w[i] << "\t";
  //std::cout  << "r: " << r << std::endl;

  for (size_t i=nsmooth; i<v.size()-nsmooth; ++i) {
    m[0] = v[i];

    //    std::cout << "before: ";
    //    for (int32_t j=0; j<nsmooth; ++j) std::cout << v[i-nsmooth+j] << "\t";
    // for (int32_t j=nsmooth; j>=0; --j) std::cout << m[j] << "\t";
    //for (int32_t j=1; j<=nsmooth; ++j) std::cout << v[i+j] << "\t";
    //std::cout << std::endl;
    double val(w[0]*m[0]);
    for (int32_t j=1; j<=nsmooth; ++j) val += w[j] * (m[j] + v[i+j]);
    v[i] = val*r;
    //std::cout << v[i] << std::endl;

    for (int32_t i=nsmooth; i>0; --i) m[i] = m[i-1];
  }
  return;
}

#endif /* _UTIL_H_ */
