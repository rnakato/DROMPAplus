/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_READFILE_H_
#define _DD_READFILE_H_

#include <unordered_map>
#include <string>
#include <vector>
#include "dd_class.hpp"
#include "SSP/common/seq.hpp"

void scan_samplestr(const std::string &str,
		    std::unordered_map<std::string, SampleFile> &sample,
		    std::vector<SamplePair> &samplepair,
		    WigType iftype);
pdSample scan_pdstr(const std::string &str);
WigArray readInputData(const std::string &filename, const int32_t binsize, const int32_t nbin, const WigType &iftype, const chrsize &chr);

#endif /* _DD_READFILE_H_ */
