#include <algorithm>
#include "alglib.h"
#include "alglib-3.10.0/src/statistics.h"
#include "alglib-3.10.0/src/specialfunctions.h"

double stdNormdist(double x, double m, double myu)
{
  double z = (x - m)/myu;
  double d = alglib::normaldistribution(z);
  double p = std::min(d,1-d)*2; // two-sided test
  return p;
}
