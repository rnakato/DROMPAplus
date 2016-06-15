/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _STATISTICS_H_
#define _STATISTICS_H_

#include <vector>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

double _getNegativeBinomial(int k, double p, double n);
double _getZIP(int k, double p, double p0);
double _getZINB(int k, double p, double n, double p0);
void iterateZINB(void *, double, double, double &, double &, double &);
double _getPoisson(int i, double m);
void iteratePoisson(void *par, double ave_pre, double &ave, double &p0);

template <class T>
T getPercentile(std::vector<T> array, double per)
{
  int binnum = array.size();
  std::vector<T> sortarray;

  for(int i=0; i<binnum; ++i) {
    if(array[i]) sortarray.push_back(array[i]);
  }
  std::sort(sortarray.begin(), sortarray.end());

  return sortarray[(int)(sortarray.size()*per)];
};

template <class T>
void getMoment(const std::vector<T> v, double &mean, double &var)
{
  mean = std::accumulate(v.begin(), v.end(), 0.0) / (double)v.size();
  var=0;
  for(auto x:v) var += (x - mean)*(x - mean);
  var /= (double)(v.size() -1);
  return;
}

#endif /* _STATISTICS_H_ */
