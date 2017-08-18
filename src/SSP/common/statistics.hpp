/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _STATISTICS_HPP_
#define _STATISTICS_HPP_

#include <vector>
#include <cmath>
#include <cstdint>

double _getNegativeBinomial(int32_t k, double p, double n);
double _getZIP(int32_t k, double p, double p0);
double _getZINB(int32_t k, double p, double n, double p0);
void iterateZINB(void *, double, double, double &, double &, double &);
double _getPoisson(int32_t i, double m);
void iteratePoisson(void *par, double ave_pre, double &ave, double &p0);
double getlogpZINB(double k, double p, double n);

namespace MyStatistics {

  template <class T>
    class moment {
    double mean;
    double var;
    double sd;
    
  public:
  moment(const std::vector<T> v, int32_t from, int32_t to=0): mean(0), var(0), sd(0) {
      if(!to) to = v.size();
      int32_t size(to - from);
      if(size) {
	for(int32_t i=from; i<to; ++i) mean += v[i];
	mean /= static_cast<double>(size);
	for(int32_t i=from; i<to; ++i) var += (v[i] - mean)*(v[i] - mean);
	var /= static_cast<double>(size -1);
	sd = sqrt(var);
      }
      //    std::cout << mean << "\t" << var << "\t" << sd << std::endl;
    }
    double getmean() const { return mean; }
    double getvar()  const { return var; }
    double getsd()   const { return sd; }
  };
  
  template <class T>
    T getPercentile(const std::vector<T> &array, const double per, int32_t binnum=0)
    {
      if(!binnum) binnum = array.size();
      std::vector<T> sortarray;
      
      for(int32_t i=0; i<binnum; ++i) {
	if(array[i]) sortarray.push_back(array[i]);
      }
      if(!sortarray.size()) return 0;
      
      sort(sortarray.begin(), sortarray.end());
      
      return sortarray[(int32_t)(sortarray.size()*per)];
    };
}

#endif /* _STATISTICS_HPP_ */
