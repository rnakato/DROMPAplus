/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */

#include <sstream>
#include <string>
#include "util.h"
#include "macro.h"

void isFile(string str)
{
  boost::filesystem::path const file(str);
  if(!boost::filesystem::exists(file)) PRINTERR(str << " does not exist.");
}

string IntToString(int n)
{
  ostringstream stream;
  stream << n;
  return stream.str();
}
