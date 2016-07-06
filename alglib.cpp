#include <algorithm>
#include "alglib.h"

double stdNormdist(double x, double m, double myu)
{
  double z = (x - m)/myu;
  double d = alglib::normaldistribution(z);
  double p = std::min(d,1-d)*2; // two-sided test
  return p;
}

double getNormdist(double s)
{
  return alglib::normaldistribution(s);
}
