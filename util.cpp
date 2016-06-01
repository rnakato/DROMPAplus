/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */

#include "util.h"
#define PRINTERR(...)   do{ cerr << "Error: " << __VA_ARGS__ << endl; exit(1); }while(0)


void isFile(string str)
{
  boost::filesystem::path const file(str);
  if(!boost::filesystem::exists(file)) PRINTERR(str << " does not exist.");
}
