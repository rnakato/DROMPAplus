/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_READFILE_H_
#define _DD_READFILE_H_

#include "dd_gv.h"

void scan_samplestr(std::string str, std::unordered_map<std::string, SampleFile> &sample, std::vector<SamplePair> &samplepair);
pdSample scan_pdstr(std::string str);
  
#endif /* _DD_READFILE_H_ */
