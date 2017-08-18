#include <iostream>
#include "statistics.h"
#include "specialfunctions.h"

int main()
{
  double p;
  double ave=10;
  for (int i=0;i<10;++i) {
    if(!i) p = alglib::poissondistribution(i, ave);
    else p = alglib::poissondistribution(i, ave) - alglib::poissondistribution(i-1, ave);
    std::cout << i << "\t" << p << std::endl;
  }
  return 0;
}
