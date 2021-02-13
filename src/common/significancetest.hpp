/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _PEAKCALL_HPP_
#define _PEAKCALL_HPP_

double getlogp_Poisson(const double val, const double myu);
double getlogp_BinomialTest(const double n1_ref, const double n2_ref, const double ratio);

#endif /* _PEAKCALL_HPP_ */
