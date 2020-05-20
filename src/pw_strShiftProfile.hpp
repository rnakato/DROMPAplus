/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _PW_STRSHIFTPROFILE_H_
#define _PW_STRSHIFTPROFILE_H_

#include <string>

class SeqStatsGenome;
class SSPstats;

int32_t setIdLongestChr(const SeqStatsGenome &genome);
void strShiftProfile(SSPstats &sspst, SeqStatsGenome &genome, const std::string &head);


#endif /*  _PW_STRSHIFTPROFILE_H_ */
