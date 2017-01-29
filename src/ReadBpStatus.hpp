/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _READBPSTATUS_HPP_
#define _READBPSTATUS_HPP_

#include <vector>
#include "SSP/common/seq.hpp"
#include "SSP/src/BpStatus.hpp"

std::vector<int32_t> readMpbl(std::string, std::string, int32_t, int32_t);
std::vector<BpStatus> readMpbl_binary(int32_t);
std::vector<BpStatus> readMpbl_binary(std::string, std::string, int32_t);
std::vector<BpStatus> arraySetBed(std::vector<BpStatus> &, std::string, const std::vector<bed> &);

#endif // _READBPSTATUS_HPP_
