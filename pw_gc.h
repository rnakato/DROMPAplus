/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_GC_H_
#define _PW_GC_H_

#include <boost/program_options.hpp>
#include "pw_gv.h"
void make_GCdist(const boost::program_options::variables_map &, Mapfile &);
void weightRead(const variables_map &values, Mapfile &p);

#endif /* _PW_GC_H_ */
