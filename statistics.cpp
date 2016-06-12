/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "statistics.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

using namespace std;

double getNegativeBinomial(int k, double p, double n)
{
  return gsl_ran_negative_binomial_pdf(k, p, n);
}

double getZINB(int k, double p, double n, double p0)
{
  double r;
  if(!k) {
    r = p0 + (1 - p0) * gsl_ran_negative_binomial_pdf(0, p, n);
  }else{
    r = (1 - p0) * gsl_ran_negative_binomial_pdf(k, p, n);
  }
  return r;
}
