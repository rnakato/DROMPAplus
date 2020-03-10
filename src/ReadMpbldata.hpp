/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _READMPBLDATA_HPP_
#define _READMPBLDATA_HPP_

#include "BpStatus.hpp"
#include "../submodules/SSP/common/BedFormat.hpp"

std::vector<int32_t> readMpblWigArray(const std::string &, const std::string &, const int32_t, const int32_t);
std::vector<BpStatus> readMpblBpArray(const std::string &, const std::string &, const int32_t, const int32_t);
void setPeak_to_MpblBpArray(std::vector<BpStatus> &array, const std::string &chrname, const std::vector<bed> &vbed);

#endif // _READMPBLDATA_HPP_
