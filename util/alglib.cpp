#include <algorithm>
#include "alglib.h"
double stdNormdist(double x, double m, double myu)
{
  double z = (x - m)/myu;
  double d = alglib::normaldistribution(z);
  double p = std::min(d,1-d)*2; // two-sided test
  return p;
}

/*double _getPoisson(int i, double m)
{
  double p;
  if(!i) p = alglib::poissondistribution(i, 1);
  else p = alglib::poissondistribution(i, 1) - alglib::poissondistribution(i-1, 1);
  return p;
}
*/
