/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_READFILE_H_
#define _DD_READFILE_H_

#include <unordered_map>
#include <string>
#include <vector>
#include "dd_class.hpp"
#include "../submodules/SSP/common/seq.hpp"

pdSample scan_pdstr(const std::string &str);

WigArray loadWigData(const std::string &filename, const SampleInfo &x, const chrsize &chr);

#endif /* _DD_READFILE_H_ */
