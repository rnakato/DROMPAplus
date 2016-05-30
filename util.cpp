/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */

#include "util.h"

void isFile(string str)
{
  boost::filesystem::path const file(str);
  if(!boost::filesystem::exists(file)) {
    cerr << "Error: " << str << " does not exist." << endl;
    exit(1);
  }
}
