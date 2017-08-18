/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _READBPSTATUS_HPP_
#define _READBPSTATUS_HPP_

#include <vector>
#include "BpStatus.hpp"

class bed;

std::vector<int32_t> readMpbl(const std::string &, const std::string &, const int32_t, const int32_t);
std::vector<BpStatus> readMpbl_binary(const std::string &, const std::string &, const int32_t);
std::vector<BpStatus> OverrideBedToArray(std::vector<BpStatus> &, const std::string &, const std::vector<bed> &);

#endif // _READBPSTATUS_HPP_
