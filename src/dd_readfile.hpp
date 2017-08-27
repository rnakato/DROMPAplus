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

void scan_samplestr(const std::string &str, const std::vector<chrsize> gt,
		    std::unordered_map<std::string, SampleFile> &sample,
		    std::vector<SamplePair> &samplepair,
		    WigType iftype);
pdSample scan_pdstr(const std::string &str);

WigArray loadWigData(const std::string &filename, const SampleFile &x, const chrsize &chr);

#endif /* _DD_READFILE_H_ */
