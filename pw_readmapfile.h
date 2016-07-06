/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _READMAPFILE_H_
#define _READMAPFILE_H_

#include <boost/program_options.hpp>
#include "seq.h"
#include "pw_gv.h"

void read_mapfile(const boost::program_options::variables_map &, Mapfile &);
void check_redundant_reads(const variables_map &values, Mapfile &p);
void estimateFragLength(const variables_map &values, Mapfile &p);

#endif /* _READMAPFILE_H_ */
