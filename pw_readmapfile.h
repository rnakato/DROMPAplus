/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _READMAPFILE_H_
#define _READMAPFILE_H_

#include <boost/program_options.hpp>
#include "util/seq.h"
#include "pw_gv.h"

void read_mapfile(const boost::program_options::variables_map &, Mapfile &, RefGenome &);

#endif /* _READMAPFILE_H_ */
