/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "statistics.hpp"

double getlogp_Poisson(const double val, const double myu)
{
  double pvalue(0);
  if (myu && val > myu) {
    pvalue = _getPoisson(val, myu);
    if (!pvalue) pvalue = 1e-300;
  }
  if (pvalue) pvalue = -log10(pvalue);
  return pvalue;
}

double getlogp_BinomialTest(const double n1_ref, const double n2_ref, const double ratio)
{
  int32_t n1, n2;
  double p(0.5);  // null model
  double pvalue(0);
  if (ratio > 1) {  /* Adjust to smaller one*/
    n1 = (int32_t)ceil(n1_ref/ratio); // rounded up
    n2 = (int32_t)ceil(n2_ref);
  } else {
    n1 = (int32_t)ceil(n1_ref);
    n2 = (int32_t)ceil(n2_ref*ratio);
  }
  if ((n1 < n2) || (!n1 && !n2)) return 0;

  pvalue = getBinomial(n1, p, n1+n2);
  if (!pvalue && n1 > n2) pvalue = 1e-300;
  if (pvalue) pvalue = -log10(pvalue);
  return pvalue;
}
