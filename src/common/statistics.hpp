/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _STATISTICS_HPP_
#define _STATISTICS_HPP_

#include <vector>
#include <cmath>
#include <cstdint>

double getNormDist(const double x, double m, double myu);
double getBinomial(int32_t k, double p, int32_t n);
double _getNegativeBinomial(int32_t k, double p, double n);
double _getZINB(int32_t k, double p, double n, double p0);
void iterateZINB(void *, double, double, double &, double &, double &);
double _getPoisson(int32_t i, double m);
void iteratePoisson(void *par, double ave_pre, double &ave, double &p0);
double getlogpZINB(double k, double p, double n);

namespace MyStatistics {

  template <class T>
  T getPercentile(const std::vector<T> &array, const double per, int32_t binnum=0)
  {
    if (!binnum) binnum = array.size();
    std::vector<T> sortarray;

    for (int32_t i=0; i<binnum; ++i) {
      if(array[i]) sortarray.push_back(array[i]);
    }
    if (!sortarray.size()) return 0;

    sort(sortarray.begin(), sortarray.end());

    return sortarray[(int32_t)(sortarray.size()*per)];
  };
}

#endif /* _STATISTICS_HPP_ */
