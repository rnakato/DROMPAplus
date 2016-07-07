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
double getlogpZINB(double k, double p, double n);

template <class T>
class moment {
  double mean;
  double var;
  double sd;

 public:
 moment(const std::vector<T> v, int from, int to=0): mean(0), var(0), sd(0) {
    if(!to) to = v.size();
    int size(to - from);
    if(size) {
      for(int i=from; i<to; ++i) mean += v[i];
      mean /= (double)size;
      for(int i=from; i<to; ++i) var += (v[i] - mean)*(v[i] - mean);
      var /= (double)(size -1);
      sd = sqrt(var);
    }
    //    std::cout << mean << "\t" << var << "\t" << sd << std::endl;
  }
  double getmean() const { return mean; }
  double getvar()  const { return var; }
  double getsd()   const { return sd; }

};

template <class T>
T getPercentile(std::vector<T> array, double per, int binnum=0)
{
  if(!binnum) binnum = array.size();
  std::vector<T> sortarray;

  for(int i=0; i<binnum; ++i) {
    if(array[i]) sortarray.push_back(array[i]);
  }
  if(!sortarray.size()) return 0;
  
  std::sort(sortarray.begin(), sortarray.end());

  return sortarray[(int)(sortarray.size()*per)];
};

#endif /* _STATISTICS_H_ */
