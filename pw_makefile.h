/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */

#ifndef _PW_MAKEFILE_H_
#define _PW_MAKEFILE_H_

#include "pw_gv.h"
#include <boost/program_options.hpp>

void makewig(const boost::program_options::variables_map &, Mapfile &);

#endif /* _PW_MAKEFILE_H_ */
