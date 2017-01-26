/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef READBPSTATUS_H
#define READBPSTATUS_H

#include <vector>
#include "SSP/src/seq.h"
#include "SSP/src/bpstatus.h"

std::vector<int32_t> readMpbl(std::string, std::string, int32_t, int32_t);
std::vector<BpStatus> readMpbl_binary(int32_t);
std::vector<BpStatus> readMpbl_binary(std::string, std::string, int32_t);
std::vector<BpStatus> arraySetBed(std::vector<BpStatus> &, std::string, const std::vector<bed> &);

#endif  // READBPSTATUS_H
