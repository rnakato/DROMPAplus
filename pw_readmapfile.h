/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _READMAPFILE_H_
#define _READMAPFILE_H_

#include "pw_gv.h"

void read_mapfile(const MyOpt::Variables &, Mapfile &);
void checkRedundantReads(const MyOpt::Variables &values, Mapfile &p);
void estimateFragLength(const MyOpt::Variables &values, Mapfile &p);

#endif /* _READMAPFILE_H_ */
