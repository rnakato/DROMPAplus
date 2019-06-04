/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */

#include "../submodules/SSP/common/statistics.hpp"

double binomial_test(const double n1_ref, const double n2_ref, const double ratio)
{
  int32_t n1, n2;
  double p(0.5);  // null model
  double pvalue;
  if (ratio > 1) {  /* 大きい方を小さい方に合わせる */
    n1 = (int)ceil(n1_ref/ratio); // rounded up
    n2 = (int)ceil(n2_ref);
  } else {
    n1 = (int)ceil(n1_ref);
    n2 = (int)ceil(n2_ref*ratio);
  }
  if((n1 < n2) || (!n1 && !n2)) return 0;
  
  pvalue = getBinomial(n1, p, n1+n2);
  if(pvalue) pvalue = -log10(pvalue);
  return pvalue;
}

