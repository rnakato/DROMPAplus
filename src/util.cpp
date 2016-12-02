/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "util.h"
#include <boost/filesystem.hpp>

std::string rmchr(const std::string &chr)
{
  std::string s;
  if(!chr.find("chr")) s = chr.substr(3);
  else s = chr;
  return s;
}

void isFile(std::string str)
{
  boost::filesystem::path const file(str);
  if(!boost::filesystem::exists(file)) {
    std::cerr << "Error: " << str << " does not exist." << std::endl;
    exit(1);
  }
}
